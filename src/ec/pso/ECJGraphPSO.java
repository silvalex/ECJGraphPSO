package ec.pso;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ec.EvolutionState;
import ec.Individual;
import ec.Problem;
import ec.simple.SimpleFitness;
import ec.simple.SimpleProblemForm;

public class ECJGraphPSO extends Problem implements SimpleProblemForm{

	@Override
	public void evaluate(final EvolutionState state, final Individual ind, final int subpopulation, final int threadnum) {
	    GraphInitializer init = (GraphInitializer) state.initializer;

		if (ind.evaluated) return;   //don't evaluate the individual if it's already evaluated
        if (!(ind instanceof GraphParticle))
            state.output.fatal("Whoa!  It's not a GraphParticle!!!",null);

        GraphParticle ind2 = (GraphParticle) ind;

        Graph graph = createNewGraph(state, init.startNode.clone(), init.endNode.clone(), init.relevant, ind2.genome);

        double a = 1.0;
        double r = 1.0;
        double t = 0.0;
        double c = 0.0;

        for (Node n : graph.nodeMap.values()) {
        	double[] qos = n.getQos();
        	a *= qos[GraphInitializer.AVAILABILITY];
        	r *= qos[GraphInitializer.RELIABILITY];
        	c += qos[GraphInitializer.COST];
        }

        // Calculate longest time
        t = findLongestPath(graph);

        a = normaliseAvailability(a, init);
        r = normaliseReliability(r, init);
        t = normaliseTime(t, init);
        c = normaliseCost(c, init);

        double fitness = init.w1 * a + init.w2 * r + init.w3 * t + init.w4 * c;

        ((SimpleFitness)ind2.fitness).setFitness(state,
                // ...the fitness...
                fitness,
                ///... is the individual ideal?  Indicate here...
                false);

        ind2.setStringRepresentation(graph.toString());
        ind2.evaluated = true;
	}

	private double normaliseAvailability(double availability, GraphInitializer init) {
		if (init.maxAvailability - init.minAvailability == 0.0)
			return 1.0;
		else
			return (availability - init.minAvailability)/(init.maxAvailability - init.minAvailability);
	}

	private double normaliseReliability(double reliability, GraphInitializer init) {
		if (init.maxReliability - init.minReliability == 0.0)
			return 1.0;
		else
			return (reliability - init.minReliability)/(init.maxReliability - init.minReliability);
	}

	private double normaliseTime(double time, GraphInitializer init) {
		if (init.maxTime - init.minTime == 0.0)
			return 1.0;
		else
			return (init.maxTime - time)/(init.maxTime - init.minTime);
	}

	private double normaliseCost(double cost, GraphInitializer init) {
		if (init.maxCost - init.minCost == 0.0)
			return 1.0;
		else
			return (init.maxCost - cost)/(init.maxCost - init.minCost);
	}

	/**
	 * Uses the Bellman-Ford algorithm with negative weights to find the longest
	 * path in an acyclic directed graph.
	 *
	 * @param g
	 * @return list of edges composing longest path
	 */
	private double findLongestPath(Graph g) {
		Map<String, Double> distance = new HashMap<String, Double>();
		Map<String, Node> predecessor = new HashMap<String, Node>();

		// Step 1: initialize graph
		for (Node node : g.nodeMap.values()) {
			if (node.getName().equals("start"))
				distance.put(node.getName(), 0.0);
			else
				distance.put(node.getName(), Double.POSITIVE_INFINITY);
		}

		// Step 2: relax edges repeatedly
		for (int i = 1; i < g.nodeMap.size(); i++) {
			for (Edge e : g.edgeList) {
				if ((distance.get(e.getFromNode().getName()) -
				        e.getToNode().getQos()[GraphInitializer.TIME])
				        < distance.get(e.getToNode().getName())) {
					distance.put(e.getToNode().getName(), (distance.get(e.getFromNode().getName()) - e.getToNode().getQos()[GraphInitializer.TIME]));
					predecessor.put(e.getToNode().getName(), e.getFromNode());
				}
			}
		}

		// Now retrieve total cost
		Node pre = predecessor.get("end");
		double totalTime = 0.0;

		while (pre != null) {
			totalTime += pre.getQos()[GraphInitializer.TIME];
			pre = predecessor.get(pre.getName());
		}

		return totalTime;
	}

	public Graph createNewGraph(EvolutionState state, Node start, Node end, Set<Node> relevant, double[] weights) {
		GraphInitializer init = (GraphInitializer) state.initializer;

		Graph newGraph = new Graph();

		Set<String> currentEndInputs = new HashSet<String>();
		Map<String,Edge> connections = new HashMap<String,Edge>();

		// Connect start node
		connectCandidateToGraphByInputs(start, connections, newGraph, currentEndInputs, init);

		Set<Node> seenNodes = new HashSet<Node>();
		List<ListItem> candidateList = new ArrayList<ListItem>();

		populateCandidateList(init.serviceToIndexMap, relevant, candidateList, weights);
		Collections.sort(candidateList);

		finishConstructingGraph(currentEndInputs, end, candidateList, connections, init, newGraph, seenNodes, relevant);

		return newGraph;
	}

	private void populateCandidateList(Map<String, Integer> serviceToIndexMap, Set<Node> relevant, List<ListItem> candidateList, double[] weights) {
		// Go through all relevant nodes
		for (Node n : relevant) {
			// Find the index for that node
			int index = serviceToIndexMap.get(n.getName());
			// Retrieve weight associated with service using the index, and create list item
			ListItem item = new ListItem(n.getName(), weights[index]);
			// Add item to list
			candidateList.add(item);
		}
	}

	private void finishConstructingGraph(Set<String> currentEndInputs, Node end, List<ListItem> candidateList, Map<String,Edge> connections,
	        GraphInitializer init, Graph newGraph, Set<Node> seenNodes, Set<Node> relevant) {

	 // While end cannot be connected to graph
		while(!checkCandidateNodeSatisfied(init, connections, newGraph, end, end.getInputs(), null)){
			connections.clear();

            // Select node
            int index;

            candidateLoop:
            for (index = 0; index < candidateList.size(); index++) {
            	ListItem item = candidateList.get(index);
            	Node candidate = init.serviceMap.get(item.serviceName).clone();
                // For all of the candidate inputs, check that there is a service already in the graph
                // that can satisfy it

                if (!checkCandidateNodeSatisfied(init, connections, newGraph, candidate, candidate.getInputs(), null)) {
                    connections.clear();
                	continue candidateLoop;
                }

                // Connect candidate to graph, adding its reachable services to the candidate list
                connectCandidateToGraphByInputs(candidate, connections, newGraph, currentEndInputs, init);
                connections.clear();

                break;
            }

            candidateList.remove(index);
        }

        connectCandidateToGraphByInputs(end, connections, newGraph, currentEndInputs, init);
        connections.clear();
        removeDanglingNodes(newGraph);
	}

	private boolean checkCandidateNodeSatisfied(GraphInitializer init,
			Map<String, Edge> connections, Graph newGraph,
			Node candidate, Set<String> candInputs, Set<Node> fromNodes) {

		Set<String> candidateInputs = new HashSet<String>(candInputs);
		Set<String> startIntersect = new HashSet<String>();

		// Check if the start node should be considered
		Node start = newGraph.nodeMap.get("start");

		if (fromNodes == null || fromNodes.contains(start)) {
    		for(String output : start.getOutputs()) {
    			Set<String> inputVals = init.taxonomyMap.get(output).servicesWithInput.get(candidate);
    			if (inputVals != null) {
    				candidateInputs.removeAll(inputVals);
    				startIntersect.addAll(inputVals);
    			}
    		}

    		if (!startIntersect.isEmpty()) {
    			Edge startEdge = new Edge(startIntersect);
    			startEdge.setFromNode(start);
    			startEdge.setToNode(candidate);
    			connections.put(start.getName(), startEdge);
    		}
		}

		for (String input : candidateInputs) {
			boolean found = false;
			for (Node s : init.taxonomyMap.get(input).servicesWithOutput) {
			    if (fromNodes == null || fromNodes.contains(s)) {
    				if (newGraph.nodeMap.containsKey(s.getName())) {
    					Set<String> intersect = new HashSet<String>();
    					intersect.add(input);

    					Edge mapEdge = connections.get(s.getName());
    					if (mapEdge == null) {
    						Edge e = new Edge(intersect);
    						e.setFromNode(newGraph.nodeMap.get(s.getName()));
    						e.setToNode(candidate);
    						connections.put(e.getFromNode().getName(), e);
    					} else
    						mapEdge.getIntersect().addAll(intersect);

    					found = true;
    					break;
    				}
			    }
			}
			// If that input cannot be satisfied, move on to another candidate
			// node to connect
			if (!found) {
				// Move on to another candidate
				return false;
			}
		}
		return true;
	}

	public void connectCandidateToGraphByInputs(Node candidate, Map<String,Edge> connections, Graph graph, Set<String> currentEndInputs, GraphInitializer init) {

		graph.nodeMap.put(candidate.getName(), candidate);
		graph.edgeList.addAll(connections.values());
		candidate.getIncomingEdgeList().addAll(connections.values());

		for (Edge e : connections.values()) {
			Node fromNode = graph.nodeMap.get(e.getFromNode().getName());
			fromNode.getOutgoingEdgeList().add(e);
		}
		for (String o : candidate.getOutputs()) {
			currentEndInputs.addAll(init.taxonomyMap.get(o).endNodeInputs);
		}
	}

	public void removeDanglingNodes(Graph graph) {
	    List<Node> dangling = new ArrayList<Node>();
	    for (Node g : graph.nodeMap.values()) {
	        if (!g.getName().equals("end") && g.getOutgoingEdgeList().isEmpty())
	            dangling.add( g );
	    }

	    for (Node d: dangling) {
	        removeDangling(d, graph);
	    }
	}

	private void removeDangling(Node n, Graph graph) {
	    if (n.getOutgoingEdgeList().isEmpty()) {
	        graph.nodeMap.remove( n.getName() );
	        for (Edge e : n.getIncomingEdgeList()) {
	            e.getFromNode().getOutgoingEdgeList().remove( e );
	            graph.edgeList.remove( e );
	            removeDangling(e.getFromNode(), graph);
	        }
	    }
	}
}
