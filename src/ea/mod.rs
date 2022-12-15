#[derive(Debug, Clone, PartialEq)]
pub struct WeightedGraph {
    edges: Vec<Vec<(usize, f64)>>,
}

impl WeightedGraph {
    pub fn len(&self) -> usize {
        self.edges.len()
    }

    pub fn new(nodes: usize) -> Self {
        Self {
            edges: vec![Vec::new(); nodes],
        }
    }

    pub fn add_edge(&mut self, edge: (usize, usize), weight: f64) {
        if let Some(edge) = self.edges[edge.0].iter_mut().find(|(v, _)| *v == edge.1) {
            edge.1 = weight;
        } else {
            self.edges[edge.0].push((edge.1, weight))
        }
    }

    pub fn weight_of(&self, edge: (usize, usize)) -> f64 {
        self.edges[edge.0]
            .iter()
            .find(|v| v.0 == edge.1)
            .unwrap_or(&(0, 0.0))
            .1
    }

    pub fn outgoing(&self, node: usize) -> &[(usize, f64)] {
        &self.edges[node]
    }

    pub fn reversed(&self) -> Self {
        let mut r = WeightedGraph::new(self.len());
        for (node, edges) in self.edges.iter().enumerate() {
            for (dest, weight) in edges {
                r.add_edge((*dest, node), *weight)
            }
        }
        r
    }
}

pub fn push_relabel(network: WeightedGraph) -> WeightedGraph {
    // define q as 0
    // define s as n-1
    let n = network.len();

    let mut f = network.clone();
    let mut ex = vec![0.0; n];
    let r = network.reversed();
    let mut h = vec![0; n];

    // create ex
    // let i = f.edges.iter().enumerate().flatten();
    // O(m)
    for (node, edges) in f.edges.iter().enumerate() {
        for (dest, weight) in edges {
            ex[node] -= weight;
            ex[*dest] += weight;
        }
    }
    ex[0] = 0.0;

    // create heights
    h[0] = n as i32;

    while let Some((u, _)) = ex.iter().skip(1).enumerate().find(|(_, v)| **v > 0.0) {
        if r.outgoing(u).iter().all(|(v, _)| h[u] <= h[*v]) {
            h[u] += 1
        } else {
            let (v, rest) = r
                .outgoing(u)
                .iter()
                .find(|(v, _)| h[*v] + 1 == h[u])
                .unwrap();

            let d = rest.min(ex[u]);
            if network.weight_of((u, *v)) != 0.0 {
                f.add_edge((u, *v), f.weight_of((u, *v)) + d)
            } else {
                f.add_edge((u, *v), f.weight_of((u, *v)) - d)
            }

            // recompute r
            todo!()
        }
    }

    f
}
