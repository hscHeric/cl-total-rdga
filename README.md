<body>
  <h1>CL-TOTAL-RDGA – Genetic Algorithm for the Total Roman Domination Problem</h1>

  <p>This repository implements a <strong>genetic algorithm</strong> for the <strong>Total Roman Domination Problem</strong> (TRDP) in C++. The solution assumes graphs without isolated vertices and uses either adjacency list (for sparse graphs) or adjacency matrix (for dense graphs) representations.</p>

  <h2>Problem Description</h2>
  <p>Given a graph G = (V, E), a function f : V → {0, 1, 2} is called a <strong>total Roman dominating function</strong> (TRDF) if:</p>
  <ul>
    <li>Every vertex v with f(v) = 0 has a neighbor u with f(u) = 2;</li>
    <li>Every vertex v with f(v) > 0 has a neighbor u with f(u) > 0.</li>
  </ul>
  <p>The goal is to minimize <code>ω(f) = ∑ f(v)</code>. The minimum value is denoted as <strong>γ<sub>tR</sub>(G)</strong>.</p>

  <h2>Dependencies</h2>
  <ul>
    <li><strong>Boost</strong>: used for <code>dynamic_bitset</code> manipulation (<code>#include &lt;boost/dynamic_bitset/dynamic_bitset.hpp&gt;</code>).</li>
    <li>A C++17-compatible compiler or newer is required.</li>
  </ul>

  <h2>Compilation</h2>
  <p>Compile the project using <code>make</code>:</p>
  <pre><code>make</code></pre>

  <p>Make sure Boost is installed. On Ubuntu, for example:</p>
  <pre><code>sudo apt install libboost-all-dev</code></pre>

  <h2>Execution</h2>
  <p>Standard mode:</p>
  <pre><code>./total_rdga &lt;path_to_graph.txt&gt; [options]</code></pre>

  <h3>Optional parameters:</h3>
  <ul>
    <li><code>--generations N</code>: maximum number of generations</li>
    <li><code>--stagnation N</code>: max generations without improvement</li>
    <li><code>--crossover VAL</code>: crossover rate</li>
    <li><code>--mutation VAL</code>: mutation rate</li>
    <li><code>--elitism VAL</code>: elitism rate</li>
    <li><code>--population VAL</code>: population factor</li>
    <li><code>--tournament N</code>: tournament size</li>
    <li><code>--trials N</code>: number of repetitions (default: 1)</li>
    <li><code>--output file.csv</code>: output CSV file path</li>
  </ul>

  <h3>Example:</h3>
  <pre><code>./total_rdga instances/graph.txt --generations 500 --trials 10</code></pre>

  <h2>Parallel Execution</h2>
  <p>Parallel execution of trials is enabled via <code>#define PARALLEL 1</code> in <code>main.cpp</code>. Default number of threads is 8:</p>
  <pre><code>#define NUM_THREADS 8</code></pre>

  <p>You can change this value directly in the source code to take advantage of multi-core machines.</p>

  <h2>Using IRACE</h2>
  <p>The algorithm supports parameter tuning using <a href="https://cran.r-project.org/web/packages/irace/index.html" target="_blank">irace</a>.</p>
  <p>To enable IRACE mode, set <code>#define IRACE 1</code> in <code>main.cpp</code> before compiling.</p>

  <h3>Execution in IRACE mode:</h3>
  <p>The program expects 4 parameters (the first 3 are ignored):</p>
  <pre><code>./total_rdga dummy dummy dummy path/graph.txt</code></pre>
  <p>Only the fitness value (a single number) will be printed to stdout.</p>

  <h3>Parameter configuration for irace:</h3>
  <pre><code>
--crossover "real (0, 1)"
--mutation "real (0, 0.3)"
--elitism "real (0, 0.5)"
--population "real (2, 10)"
--tournament "integer (2, 10)"
--generations "integer (100, 1000)"
--stagnation "integer (100, 1000)"
</code></pre>

  <h2>Input Example</h2>
  <p>The graph must be in a <code>.txt</code> file in edge list format (no weights):</p>
  <pre><code>
0 1
0 2
1 3
2 3
  </code></pre>

  <h2>Heuristics</h2>
  <p>The initial population is generated using 5 distinct heuristics that guarantee valid starting solutions, as described in the submitted paper.</p>

  <h2>License</h2>
  <p>This project is licensed under the <strong>GNU General Public License v2.0</strong>.</p>
  <p>You are free to redistribute and/or modify it under the terms published by the Free Software Foundation.</p>
  <p>Full license text: <a href="https://www.gnu.org/licenses/old-licenses/gpl-2.0.html" target="_blank">https://www.gnu.org/licenses/old-licenses/gpl-2.0.html</a></p>

  <hr>
  <p>Heric da Silva Cruz – <em>hericsilvaho@gmail.com</em></p>
  <p>Prof. Atílio Gomes Luiz – <em>gomes.atilio@ufc.br</em></p>
</body>
