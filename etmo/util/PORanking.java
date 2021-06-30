//  Ranking.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package etmo.util;

import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.util.comparators.OverallConstraintViolationComparator;
import etmo.util.comparators.PODominanceComparator;

import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * This class implements some facilities for ranking solutions. Given a
 * <code>SolutionSet</code> object, their solutions are ranked according to
 * scheme proposed in NSGA-II; as a result, a set of subsets are obtained. The
 * subsets are numbered starting from 0 (in NSGA-II, the numbering starts from
 * 1); thus, subset 0 contains the non-dominated solutions, subset 1 contains
 * the non-dominated solutions after removing those belonging to subset 0, and
 * so on.
 */
public class PORanking {

	/**
	 * The <code>SolutionSet</code> to rank
	 */
	private SolutionSet solutionSet_;

	/**
	 * An array containing all the fronts found during the search
	 */
	private SolutionSet[] ranking_;

	/**
	 * stores a <code>Comparator</code> for dominance checking
	 */
	private static Comparator dominance_;

	/**
	 * stores a <code>Comparator</code> for Overal Constraint Violation
	 * Comparator checking
	 */
	private static Comparator constraint_ = new OverallConstraintViolationComparator();

	public PORanking(SolutionSet solutionSet, boolean[] isChosen) {
		solutionSet_ = solutionSet;

		dominance_ = new PODominanceComparator(isChosen);

		// dominateMe[i] contains the number of solutions dominating i
		int[] dominateMe = new int[solutionSet_.size()];

		// iDominate[k] contains the list of solutions dominated by k
		List<Integer>[] iDominate = new List[solutionSet_.size()];

		// front[i] contains the list of individuals belonging to the front i
		List<Integer>[] front = new List[solutionSet_.size() + 1];

		// flagDominate is an auxiliar encodings.variable
		int flagDominate;

		// Initialize the fronts
		for (int i = 0; i < front.length; i++)
			front[i] = new LinkedList<Integer>();

		/*
		 * //-> Fast non dominated sorting algorithm for (int p = 0; p <
		 * solutionSet_.size(); p++) { // Initialice the list of individuals
		 * that i dominate and the number // of individuals that dominate me
		 * iDominate[p] = new LinkedList<Integer>(); dominateMe[p] = 0; // For
		 * all q individuals , calculate if p dominates q or vice versa for (int
		 * q = 0; q < solutionSet_.size(); q++) { flagDominate
		 * =constraint_.compare(solutionSet.get(p),solutionSet.get(q)); if
		 * (flagDominate == 0) { flagDominate
		 * =dominance_.compare(solutionSet.get(p),solutionSet.get(q)); }
		 * 
		 * if (flagDominate == -1) { iDominate[p].add(new Integer(q)); } else if
		 * (flagDominate == 1) { dominateMe[p]++; } }
		 * 
		 * // If nobody dominates p, p belongs to the first front if
		 * (dominateMe[p] == 0) { front[0].add(new Integer(p));
		 * solutionSet.get(p).setRank(0); } }
		 */

		// -> Fast non dominated sorting algorithm
		// Contribution of Guillaume Jacquenot
		for (int p = 0; p < solutionSet_.size(); p++) {
			// Initialize the list of individuals that i dominate and the number
			// of individuals that dominate me
			iDominate[p] = new LinkedList<Integer>();
			dominateMe[p] = 0;
		}
		for (int p = 0; p < (solutionSet_.size() - 1); p++) {
			// For all q individuals , calculate if p dominates q or vice versa
			for (int q = p + 1; q < solutionSet_.size(); q++) {
				flagDominate = constraint_.compare(solutionSet.get(p), solutionSet.get(q));
				if (flagDominate == 0) {
					flagDominate = dominance_.compare(solutionSet.get(p), solutionSet.get(q));
				}
				if (flagDominate == -1) {
					iDominate[p].add(q);
					dominateMe[q]++;
				} else if (flagDominate == 1) {
					iDominate[q].add(p);
					dominateMe[p]++;
				}
			}
			// If nobody dominates p, p belongs to the first front
		}
		for (int p = 0; p < solutionSet_.size(); p++) {
			if (dominateMe[p] == 0) {
				front[0].add(p);
				solutionSet.get(p).setRank(0);
			}
		}

		// Obtain the rest of fronts
		int i = 0;
		Iterator<Integer> it1, it2; // Iterators
		while (front[i].size() != 0) {
			i++;
			it1 = front[i - 1].iterator();
			while (it1.hasNext()) {
				it2 = iDominate[it1.next()].iterator();
				while (it2.hasNext()) {
					int index = it2.next();
					dominateMe[index]--;
					if (dominateMe[index] == 0) {
						front[i].add(index);
						solutionSet_.get(index).setRank(i);
					}
				}
			}
		}
		// <-
		//将非支配排序的个体按照索引存如ranking_中
		ranking_ = new SolutionSet[i];
		// 0,1,2,....,i-1 are front, then i fronts
		for (int j = 0; j < i; j++) {
			ranking_[j] = new SolutionSet(front[j].size());
			it1 = front[j].iterator();
			while (it1.hasNext()) {
				ranking_[j].add(solutionSet.get(it1.next()));
			}
		}

	}

	/**
	 * Constructor.
	 * 
	 * @param solutionSet
	 *            The <code>SolutionSet</code> to be ranked.
	 */
	public PORanking(SolutionSet solutionSet, Solution objSet) {
		solutionSet_ = solutionSet;

		dominance_ = new PODominanceComparator(objSet);

		// dominateMe[i] contains the number of solutions dominating i
		int[] dominateMe = new int[solutionSet_.size()];

		// iDominate[k] contains the list of solutions dominated by k
		List<Integer>[] iDominate = new List[solutionSet_.size()];

		// front[i] contains the list of individuals belonging to the front i
		List<Integer>[] front = new List[solutionSet_.size() + 1];

		// flagDominate is an auxiliar encodings.variable
		int flagDominate;

		// Initialize the fronts
		for (int i = 0; i < front.length; i++)
			front[i] = new LinkedList<Integer>();

		/*
		 * //-> Fast non dominated sorting algorithm for (int p = 0; p <
		 * solutionSet_.size(); p++) { // Initialice the list of individuals
		 * that i dominate and the number // of individuals that dominate me
		 * iDominate[p] = new LinkedList<Integer>(); dominateMe[p] = 0; // For
		 * all q individuals , calculate if p dominates q or vice versa for (int
		 * q = 0; q < solutionSet_.size(); q++) { flagDominate
		 * =constraint_.compare(solutionSet.get(p),solutionSet.get(q)); if
		 * (flagDominate == 0) { flagDominate
		 * =dominance_.compare(solutionSet.get(p),solutionSet.get(q)); }
		 * 
		 * if (flagDominate == -1) { iDominate[p].add(new Integer(q)); } else if
		 * (flagDominate == 1) { dominateMe[p]++; } }
		 * 
		 * // If nobody dominates p, p belongs to the first front if
		 * (dominateMe[p] == 0) { front[0].add(new Integer(p));
		 * solutionSet.get(p).setRank(0); } }
		 */

		// -> Fast non dominated sorting algorithm
		// Contribution of Guillaume Jacquenot
		for (int p = 0; p < solutionSet_.size(); p++) {
			// Initialize the list of individuals that i dominate and the number
			// of individuals that dominate me
			iDominate[p] = new LinkedList<Integer>();
			dominateMe[p] = 0;
		}
		for (int p = 0; p < (solutionSet_.size() - 1); p++) {
			// For all q individuals , calculate if p dominates q or vice versa
			for (int q = p + 1; q < solutionSet_.size(); q++) {
				flagDominate = constraint_.compare(solutionSet.get(p), solutionSet.get(q));
				if (flagDominate == 0) {
					flagDominate = dominance_.compare(solutionSet.get(p), solutionSet.get(q));
				}
				if (flagDominate == -1) {
					iDominate[p].add(q);
					dominateMe[q]++;
				} else if (flagDominate == 1) {
					iDominate[q].add(p);
					dominateMe[p]++;
				}
			}
			// If nobody dominates p, p belongs to the first front
		}
		for (int p = 0; p < solutionSet_.size(); p++) {
			if (dominateMe[p] == 0) {
				front[0].add(p);
				solutionSet.get(p).setRank(0);
			}
		}
		/*按第0层举例,当i=0时，判断第0层的非支配个体数，若不为0，则取迭代器，it1,若it1.hasNext()不为0，则it1.next()为第0层的第一个个体下标索引，
		由于每个个体都用list存放该个体所支配的个体，则取迭代器it2,若it2.hasNext()不为0，则it2.next()为被it1.next()个体支配的个体的下标索引，
		dominateME[index]--表示支配index的个体数减1，若上一层非支配个体循环完后，如果dominateME[index]为0，则代表该index索引的个体属于下一层
		非支配个体
		* */
		// Obtain the rest of fronts
		int i = 0;
		Iterator<Integer> it1, it2; // Iterators
		while (front[i].size() != 0) {
			i++;
			it1 = front[i - 1].iterator();
			while (it1.hasNext()) {
				it2 = iDominate[it1.next()].iterator();
				while (it2.hasNext()) {
					int index = it2.next();
					dominateMe[index]--;
					if (dominateMe[index] == 0) {
						front[i].add(index);
						solutionSet_.get(index).setRank(i);
					}
				}
			}
		}
		// <-

		ranking_ = new SolutionSet[i];
		// 0,1,2,....,i-1 are front, then i fronts
		for (int j = 0; j < i; j++) {
			ranking_[j] = new SolutionSet(front[j].size());
			it1 = front[j].iterator();
			while (it1.hasNext()) {
				ranking_[j].add(solutionSet.get(it1.next()));
			}
		}

	} // Ranking

	/**
	 * Returns a <code>SolutionSet</code> containing the solutions of a given
	 * rank.
	 * 
	 * @param rank
	 *            The rank
	 * @return Object representing the <code>SolutionSet</code>.
	 */
	public SolutionSet getSubfront(int rank) {
		return ranking_[rank];
	} // getSubFront

	/**
	 * Returns the total number of subFronts founds.
	 */
	public int getNumberOfSubfronts() {
		return ranking_.length;
	} // getNumberOfSubfronts
} // Ranking