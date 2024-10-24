#![feature(generic_const_exprs)]

use std::{arch::aarch64::int16x4_t, collections::{HashMap, HashSet}};

use matrix::Matrix;
use matrix_kit::*;

use minilp::{ComparisonOp, OptimizationDirection, Problem, Variable};

fn bp_to_set<const N: usize>(bits: usize) -> Vec<usize> {
    let mut set = Vec::<usize>::new();
    for v in 0..N {
        // check if v is in this vertex selection
        if bits & (1 << v) != 0 {
            set.push(v);
        }
    }
    set
}

fn set_to_bp(set: Vec<usize>) -> usize {
    let mut bit_pattern = 0;
    for v in set {
        bit_pattern |= 1 << v;
    }
    bit_pattern
}

/// An undirected graph with N nodes
#[derive(Clone, Copy)]
struct Graph<const N: usize> where [() ; N * N]: Sized {
    pub adjacencies: Matrix<N, N, i8>
}

impl<const N: usize> Graph<N> where [(); N * N]: {

    /// Creates a graph from a list of adjacencies
    fn from_adjacency_list(adj_list: Vec<(usize, usize)>) -> Graph<N> {
        let mut adj_mat = Matrix::<N, N, i8>::new();
        for (u, v) in adj_list {
            adj_mat[u][v] = 1;
            adj_mat[v][u] = 1;
        }
        Graph { adjacencies: adj_mat }
    }

    /// Computes the complement of this graph
    fn complement(self) -> Graph<N> {
        let all_ones = Matrix::<N, N, i8>::from_flatmap([1 ; N * N]);
        let complememt_adj_mat = all_ones - self.adjacencies;
        Graph { adjacencies: complememt_adj_mat }
    }

    /// Inefficiently computes all cliques in a graph,
    /// returning a set of sets
    fn cliques(self) -> Vec<Vec<usize>> {
        // Iterate through every possible subset of vertices, and check if its a clique. If it is,
        // add that subset to our list!
        let mut cliques = Vec::<Vec<usize>>::new();

        for subset_bitmap in 1..(1 << N) { // we start at 1 because we don't care about 0 which represents choosing no vertices
            let selected_vertices = bp_to_set::<N>(subset_bitmap);

            // Now, check that every selected vertex is adjacent to every other selected vertex.
            let mut is_clique = true;
            
            for v in selected_vertices.clone() { // Check vertex v
                for u in selected_vertices.clone() { 
                    if v == u { continue; } // We don't care about self-adjacency

                    // Check that v is adjacent to u
                    if self.adjacencies[v][u] == 0 {
                        // this is not a clique!
                        is_clique = false;
                    }
                }
            }

            if is_clique {
                cliques.push(selected_vertices);
            }
        }

        cliques
    }

    /// Finds all independent sets in a graph
    fn independent_sets(self) -> Vec<Vec<usize>> {
        self.complement().cliques()
    }

    /// Finds all independent sets that contain a specific vertex
    fn independent_sets_containing(self, v: usize) -> Vec<Vec<usize>> {
        self.independent_sets()
            .into_iter()
            .filter(|is| is.contains(&v))
            .collect()
    }
}

fn main() {
    let graph = Graph::<6>::from_adjacency_list(vec![
        (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
        (1, 2), (2, 3), (3, 4), (4, 5), (5, 1)
    ]);

    let mut problem = Problem::new(OptimizationDirection::Minimize);

    // we will have a "vector" of variables, each associated with an independent set of vertices.
    let mut vars = HashMap::<usize, Variable>::new();
    
    for independent_set in graph.independent_sets() {
        // We will turn this independent set into a bit pattern where a 1 means we choose that vertex
        let bit_pattern = set_to_bp(independent_set);
        vars.insert(bit_pattern, problem.add_var(1.0, (0.0, f64::INFINITY))); // Our objective function is just the sum of our variables
    }
    
    // Now, we add the constraint to make sure every vertex is covered.
    for v in 0..6 {
        let ind_sets_with_v = graph.independent_sets_containing(v);
        println!("Adding sum for node {:?}", v);

        let mut sum = Vec::new();

        for ind_set in ind_sets_with_v {
            println!("\tAdding term to sum for node {:?}: x_{:?}", v, ind_set);
            sum.push((*vars.get(&set_to_bp(ind_set)).unwrap(), 1.0));
            println!("\t\tPushed: {:?}", sum.last())
        }

        problem.add_constraint(sum, ComparisonOp::Ge, 1.0);
    }

    let solution = problem.solve().unwrap();
    println!("Optimal Value: {:?}", solution.objective());
}