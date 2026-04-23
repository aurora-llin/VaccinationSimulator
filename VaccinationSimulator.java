import java.util.*;

/**
 * Vaccination Strategy Simulation
 *
 * Compare:
 * 3 Network Types:
 *   1. Scale-Free (Barabasi-Albert)
 *   2. Random Graph (Erdos-Renyi)
 *   3. Small World (Watts-Strogatz)
 *
 * 3 Vaccination Strategies:
 *   1. Degree Centrality
 *   2. Betweenness Centrality
 *   3. Random
 *
 * Vaccination Rates:
 *   5%,15%,30%, 50%
 */

public class VaccinationSimulator {

    // =====================================================
    // GRAPH
    // =====================================================
    static class Graph {
        int n;
        List<List<Integer>> adj;

        Graph(int n) {
            this.n = n;
            adj = new ArrayList<>();
            for (int i = 0; i < n; i++)
                adj.add(new ArrayList<>());
        }

        void addEdge(int u, int v) {
            if (u == v) return;
            if (!adj.get(u).contains(v)) {
                adj.get(u).add(v);
                adj.get(v).add(u);
            }
        }

        void removeEdge(int u, int v) {
            adj.get(u).remove(Integer.valueOf(v));
            adj.get(v).remove(Integer.valueOf(u));
        }

        int degree(int u) {
            return adj.get(u).size();
        }

        int maxDegree() {
            int max = 0;
            for (int i = 0; i < n; i++)
                max = Math.max(max, degree(i));
            return max;
        }

        double avgDegree() {
            long total = 0;
            for (int i = 0; i < n; i++)
                total += degree(i);
            return (double) total / n;
        }
    }

    // =====================================================
    // NETWORK MODELS
    // =====================================================

    // 1. Scale-Free
    static Graph barabasiAlbert(int n, int m, Random rng) {
        Graph g = new Graph(n);
        int m0 = m + 1;

        for (int i = 0; i < m0; i++)
            for (int j = i + 1; j < m0; j++)
                g.addEdge(i, j);

        for (int newNode = m0; newNode < n; newNode++) {
            List<Integer> pool = new ArrayList<>();

            for (int i = 0; i < newNode; i++)
                for (int d = 0; d < g.degree(i); d++)
                    pool.add(i);

            Set<Integer> chosen = new HashSet<>();

            while (chosen.size() < m && !pool.isEmpty())
                chosen.add(pool.get(rng.nextInt(pool.size())));

            for (int target : chosen)
                g.addEdge(newNode, target);
        }

        return g;
    }

    // 2. Random Graph
    static Graph erdosRenyi(int n, double p, Random rng) {
        Graph g = new Graph(n);

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (rng.nextDouble() < p)
                    g.addEdge(i, j);
            }
        }

        return g;
    }

    // 3. Small World
    static Graph wattsStrogatz(int n, int k, double beta, Random rng) {
        Graph g = new Graph(n);

        // Ring lattice
        for (int i = 0; i < n; i++) {
            for (int j = 1; j <= k / 2; j++) {
                g.addEdge(i, (i + j) % n);
            }
        }

        // Rewire
        for (int i = 0; i < n; i++) {
            for (int j = 1; j <= k / 2; j++) {
                int neighbor = (i + j) % n;

                if (rng.nextDouble() < beta) {
                    g.removeEdge(i, neighbor);

                    int newNeighbor;
                    do {
                        newNeighbor = rng.nextInt(n);
                    } while (newNeighbor == i || g.adj.get(i).contains(newNeighbor));

                    g.addEdge(i, newNeighbor);
                }
            }
        }

        return g;
    }

    // =====================================================
    // VACCINATION STRATEGIES
    // =====================================================

    static Set<Integer> degreeVaccination(Graph g, int k) {
        PriorityQueue<int[]> pq =
                new PriorityQueue<>((a, b) -> b[1] - a[1]);

        for (int i = 0; i < g.n; i++)
            pq.offer(new int[]{i, g.degree(i)});

        Set<Integer> result = new HashSet<>();

        while (result.size() < k)
            result.add(pq.poll()[0]);

        return result;
    }

    static Set<Integer> randomVaccination(Graph g, int k, Random rng) {
        List<Integer> nodes = new ArrayList<>();

        for (int i = 0; i < g.n; i++)
            nodes.add(i);

        Collections.shuffle(nodes, rng);

        return new HashSet<>(nodes.subList(0, k));
    }

    static Set<Integer> betweennessVaccination(Graph g, int k) {
        double[] score = new double[g.n];

        for (int s = 0; s < g.n; s++) {

            Deque<Integer> stack = new ArrayDeque<>();
            List<List<Integer>> pred = new ArrayList<>();

            for (int i = 0; i < g.n; i++)
                pred.add(new ArrayList<>());

            int[] dist = new int[g.n];
            Arrays.fill(dist, -1);

            double[] sigma = new double[g.n];
            sigma[s] = 1;
            dist[s] = 0;

            Queue<Integer> q = new LinkedList<>();
            q.add(s);

            while (!q.isEmpty()) {
                int v = q.poll();
                stack.push(v);

                for (int w : g.adj.get(v)) {
                    if (dist[w] < 0) {
                        q.add(w);
                        dist[w] = dist[v] + 1;
                    }

                    if (dist[w] == dist[v] + 1) {
                        sigma[w] += sigma[v];
                        pred.get(w).add(v);
                    }
                }
            }

            double[] delta = new double[g.n];

            while (!stack.isEmpty()) {
                int w = stack.pop();

                for (int v : pred.get(w))
                    delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w]);

                if (w != s)
                    score[w] += delta[w];
            }
        }

        Integer[] nodes = new Integer[g.n];
        for (int i = 0; i < g.n; i++) nodes[i] = i;

        Arrays.sort(nodes, (a, b) -> Double.compare(score[b], score[a]));

        Set<Integer> result = new HashSet<>();

        for (int i = 0; i < k; i++)
            result.add(nodes[i]);

        return result;
    }

    // =====================================================
    // DISEASE SPREAD BFS
    // =====================================================

    static class Result {
        int infected;
        double percent;
    }

    static Result simulate(Graph g, Set<Integer> vaccinated, Random rng) {

        int[] state = new int[g.n];
        for (int v : vaccinated)
            state[v] = 3;

        int seed;
        do {
            seed = rng.nextInt(g.n);
        } while (state[seed] == 3);

        Queue<Integer> q = new LinkedList<>();
        q.add(seed);
        state[seed] = 1;

        int infected = 1;

        while (!q.isEmpty()) {
            int size = q.size();

            for (int i = 0; i < size; i++) {
                int curr = q.poll();

                for (int next : g.adj.get(curr)) {
                    if (state[next] == 0) {
                        state[next] = 1;
                        infected++;
                        q.add(next);
                    }
                }

                state[curr] = 2;
            }
        }

        Result r = new Result();
        r.infected = infected;
        r.percent = 100.0 * infected / g.n;

        return r;
    }

    // =====================================================
    // MAIN
    // =====================================================

    public static void main(String[] args) {

        Scanner sc = new Scanner(System.in);
        Random rng = new Random(42);

        boolean running = true;

        while (running) {

            System.out.println("\n======================================");
            System.out.println(" Vaccination Strategy Simulation ");
            System.out.println("======================================");
            System.out.println("1. Small Example");
            System.out.println("2. Large Example");
            System.out.println("3. Custom Input");
            System.out.println("4. Exit");
            System.out.print("Enter choice: ");

            int choice = sc.nextInt();

            if (choice == 4) {
                System.out.println("Exiting simulation...");
                break;
            }

            int N;
            String networkType;
            double vaccinationRate;
            String strategy;

            switch (choice) {
                case 1:
                    N = 100;
                    networkType = "Scale-Free";
                    vaccinationRate = 0.10;
                    strategy = "Degree";
                    break;

                case 2:
                    N = 1000;
                    networkType = "Scale-Free";
                    vaccinationRate = 0.15;
                    strategy = "Degree";
                    break;

                case 3:
                    System.out.print("Enter nodes: ");
                    N = sc.nextInt();

                    sc.nextLine(); // clear newline

                    System.out.print("Enter network type (Scale-Free / Random / Small-World): ");
                    networkType = sc.nextLine();

                    System.out.print("Enter vaccination rate (%): ");
                    vaccinationRate = sc.nextDouble() / 100.0;

                    sc.nextLine(); // clear newline

                    System.out.print("Enter strategy (Degree / Betweenness / Random): ");
                    strategy = sc.nextLine();
                    break;

                default:
                    System.out.println("Invalid choice. Try again.");
                    continue;
            }

            Graph g;

            if (networkType.equalsIgnoreCase("Scale-Free")) {
                g = barabasiAlbert(N, 2, rng);
            } else if (networkType.equalsIgnoreCase("Random")) {
                g = erdosRenyi(N, 0.03, rng);
            } else {
                g = wattsStrogatz(N, 6, 0.20, rng);
            }

            int k = (int) Math.ceil(vaccinationRate * N);

            Set<Integer> vaccinated;

            if (strategy.equalsIgnoreCase("Degree")) {
                vaccinated = degreeVaccination(g, k);
            } else if (strategy.equalsIgnoreCase("Betweenness")) {
                vaccinated = betweennessVaccination(g, k);
            } else {
                vaccinated = randomVaccination(g, k, rng);
            }
            //average result after 30 runs
            int runs = 30;
            double total = 0;
            Result result = null;

            for (int i = 0; i < runs; i++) {
                result = simulate(g, vaccinated, rng);
                total += result.percent;
            }

            double avg = total / runs;
            result.percent = avg;
            result.infected = (int) Math.round(avg * N / 100.0);

            System.out.println("Average infected: " + avg + "%");

            System.out.println("\n======================================");
            System.out.println("Simulation Result");
            System.out.println("======================================");
            System.out.println("Nodes: " + N);
            System.out.println("Network: " + networkType);
            System.out.println("Vaccination Rate: " + (vaccinationRate * 100) + "%");
            System.out.println("Strategy: " + strategy);

            printResult(strategy, result, N);

            System.out.println("\nRun another simulation...");
        }

        sc.close();
    }

    static void printResult(String name, Result r, int N) {
        System.out.printf("%-14s Infected: %3d/%d (%.1f%%)%s%n",
                name,
                r.infected,
                N,
                r.percent,
                r.percent < 20 ? "  <-- CONTAINED" : "");
    }
}