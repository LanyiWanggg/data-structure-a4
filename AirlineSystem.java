import java.io.*;
import java.util.*;
import java.util.Iterator;
import java.util.NoSuchElementException;


public class AirlineSystem
{
    private static String[] cityNames = null;                      //string storing cities
    private static Digraph G = null;                               //graph
    public static Scanner scan = null;
    private static final int INFINITY = Integer.MAX_VALUE;
    private static String fileName = "";
    static AirlineSystem airline = new AirlineSystem();
    private static int shortestDistanceViaThirdCity = 0;
    
    public static void main(String[] args) throws Exception
    {
        scan = new Scanner(System.in);
        airline.readGraph();                                //read data input file into a graph
        while(true)
        {
            switch(airline.menu())                          //show user the menu options
            {
                case 1:
                    airline.printGraph();                   //display all routes
                    break;
                case 2:
                    airline.displayMinimumSpanningTree();   //display a minimum spanning tree
                    break;
                case 3:
                    airline.getShortestPathByDistance();    
                    break;
                case 4:
                    airline.getShortestPathByCost();
                    break;
                case 5:
                    airline.getShortestPathByHops();
                    break;
                case 6:
                    airline.getAllTripsUnder();             //print out all paths under a specific amount
                    break;
                case 7:
                    airline.addRoute();                     //add a route between 2 existed cities
                    break;
                case 8:
                    airline.removeRoute();                  //remove a route between 2 existed cities
                    break;
                case 9:
                    airline.quitAndSave();                  //rewrite the graph into input file
                    break;
                case 10:
                    airline.findShortestPathThruACity();    //extra credit
                    break;
                default:
                    System.out.println("Invalid option. Please enter a correct option(1-10): ");

            }
        }
    }

    private int menu()
    {
        System.out.println("********************************* MAIN MENU *********************************");
        System.out.println("1. Show the entire list of direct routes, distances and prices.");
        System.out.println("2. Display a minimum spanning tree for the service routes based on distances.");
        System.out.println("3. Find the trip with shortest distance between two cities.");
        System.out.println("4. Find the trip with lowest cost between two cities.");
        System.out.println("5. Find the trip with least number of hops between two cities.");
        System.out.println("6. Print out all trips whose cost is less than or equal to a certain cost.");
        System.out.println("7. Add a new route to schedule.");
        System.out.println("8. Remove a route from schedule.");
        System.out.println("9. Save all routes back to the file and quit the program.");
        System.out.println("10. Find the shortest path from source city to destination city thru a third city.");
        System.out.println("*****************************************************************************");
        System.out.print("Please choose a menu option (1-10): ");

        int choice = Integer.parseInt(scan.nextLine());
        return choice;
    }

    private void readGraph() throws IOException
    {
        System.out.println("Please enter graph file name: ");
        fileName = scan.nextLine();
        Scanner fileScan = new Scanner(new FileInputStream(fileName));
        int v = fileScan.nextInt();                         //store number of cities into v
        G = new Digraph(v);                                 //create graph with v vertices

        cityNames = new String[v];                          //initialize the array storing all cities
        String st = fileScan.nextLine();
        for(int i = 0; i<v; i++)                            //store v cities into cityNames string
        {
            cityNames[i] = fileScan.nextLine();
        }

        while(fileScan.hasNext())                           //initialize graph: put cities into vertices and build edges 
        {                                                   //with cost and distance as weights of edge
            int from = fileScan.nextInt();
            int to = fileScan.nextInt();
            int distance = fileScan.nextInt();
            double cost = fileScan.nextDouble();
            G.addEdge(new DirectedEdge(from-1, to-1, cost, distance));
            G.addEdge(new DirectedEdge(to-1, from-1, cost, distance));
            //fileScan.nextLine();
        }
        fileScan.close();
        System.out.println("Data imported successfully.");
    }
    
    private void printGraph()
    {
        for(int i = 0; i<G.v; i++)                          //iterate thru vertices of graph
        {
            System.out.println(cityNames[i] + " to: ");     //print out each vertice(city)
            for(DirectedEdge e : G.adj(i))                  //iterate thru edges(cost and distance) of each vetice
            {
                System.out.println("    " + cityNames[e.to()] + " (" + e.distance() + " miles) $" + e.cost());
            }
            System.out.println();
        }
    }

    private void displayMinimumSpanningTree()
    {
        PrimMSTTrace mst = new PrimMSTTrace(G);             //create a mimimum spanning tree and initialize it  
                                                            //by putting graph into
        int v, w;                                           //v:start city  w:end city
        for(DirectedEdge e : mst.edges())                   //iterate edges of mst
        {
            v = e.either();                                 //either() and other() are two sides of an edge
            w = e.other(v);
            System.out.println(cityNames[v] + "," + cityNames[w] + ": " + e.distance());
        }
        System.out.println();
    }

    private void getShortestPathByDistance()                //from lab10
    {
        System.out.print("Please enter source city (1-" + cityNames.length + "): ");
        int source = Integer.parseInt(scan.nextLine());
        System.out.print("Please enter destination city (1-" + cityNames.length + "): ");
        int destination = Integer.parseInt(scan.nextLine());
       
        source--;                                         //decrement both since cities starts from 0 in code
        destination--;                                    //but starts from 1 from users' point of view

        G.digkstras(source, destination, false);          //the 3rd parameter is a boolean; false means by distance
                                                          //true means by cost
        if(!G.marked[destination])                        //marked array check is there is an edge between
        {                                                 //source and destination
            System.out.println("There is no route from " + cityNames[source] + " to " + cityNames[destination]);
        }
        else
        {
            Stack<Integer> path = new Stack<>();          //build a stack to store all the possible cities between s and destination
            for(int x = destination; x!=source; x=G.edgeTo[x]) //go to a previous edge on shortest path each time
            {
                path.push(x);
            }

            System.out.println("SHORTEST DISTANCE PATH from "+cityNames[source] + " to " + cityNames[destination] + " is " + G.distTo[destination] + " miles");

            int preVertex = source;                       //pop those cities and calculate the distance and print out
            System.out.print(cityNames[source] + " ");
            while (!path.empty()) {
                int v = path.pop();
                System.out.print(G.distTo[v] - G.distTo[preVertex] + " " + cityNames[v] + " ");
                preVertex = v;
              }
              System.out.println();
        }
        System.out.println();

    }

    private void getShortestPathByCost()                    //from lab10
    {                                                       //this is very similar to the distance one, just replace several words
                                                            //since we deal with the main difference between cost and distance in 
                                                            //digkstras() method in Digraph class
        System.out.print("Please enter source city (1-" + cityNames.length + "): ");
        int source = Integer.parseInt(scan.nextLine());
        System.out.print("Please enter destination city (1-" + cityNames.length + "): ");
        int destination = Integer.parseInt(scan.nextLine());

        source--;
        destination--;

        G.digkstras(source, destination, true);

        if(!G.marked[destination])
        {
            System.out.println("There is no route from " + cityNames[source] + " to " + cityNames[destination]);
        }
        else
        {
            Stack<Integer> path = new Stack<>();
            for(int x = destination; x!=source; x=G.edgeTo[x])
            {
                path.push(x);
            }

            System.out.println("LOWESR COST PATH from "+cityNames[source] + " to " + cityNames[destination] + " is $" + G.distTo[destination]);

            int preVertex = source;
            System.out.print(cityNames[source] + " ");
            while (!path.empty()) {
                int v = path.pop();
                System.out.print(G.distTo[v] - G.distTo[preVertex] + " " + cityNames[v] + " ");
                preVertex = v;
              }
              System.out.println();
        }
        System.out.println();
    }

    private void getShortestPathByHops()                    //from lab10
    {                                                       //this is a bit ddifferent than the above two since its using bfs
        System.out.print("Please enter source city (1-" + cityNames.length + "): ");
        int source = Integer.parseInt(scan.nextLine());
        System.out.print("Please enter destination city (1-" + cityNames.length + "): ");
        int destination = Integer.parseInt(scan.nextLine());

        source--;
        destination--;

        G.bfs(source);

        if(!G.marked[destination])
        {
            System.out.println("There is no route from " + cityNames[source] + " to " + cityNames[destination]);
        }
        else
        {
            Stack<Integer> path = new Stack<>();
            for(int x = destination; x!=source; x=G.edgeTo[x])
            {
                path.push(x);
            }
            path.push(source);

            System.out.println("FEWEST hops from "+cityNames[source] + " to " + cityNames[destination] + " is " + G.distTo[destination]);
 
            while (!path.empty()) 
            {
                System.out.print(cityNames[path.pop()] + " ");
            }
            System.out.println();
        }
        System.out.println();
    }

    private void getAllTripsUnder()
    {
        System.out.print("Please enter a maximum cost for any flight: ");
        double maxCost = Double.parseDouble(scan.nextLine());     
        System.out.println("ALL PATHS OF COST " + maxCost+  " OR LESS");
        System.out.println("Note that paths are duplicated, once from each end city's point of view");
        System.out.println("List of paths at most " + maxCost+  " in length: ");
        System.out.println("-----------------------------------------------");
        // for each city
        for(int i = 0; i<G.v; i++)
        {
            // run DepthFirstAmountPaths from the city and specify the amount
            DepthFirstAmountPaths dfa = new DepthFirstAmountPaths(G, i, maxCost);
            // store the results of the DFS into an array of linkedlists
            ArrayList<LinkedList<DirectedEdge>> alist = dfa.getResults();
            // iterate over each linkedlist of results in the array
            for(LinkedList<DirectedEdge> linkedList : alist)
            {
                // initialize a tripCost variable
                double tripCost=0;
                // iterate over each edge in the linkedlist
                for(DirectedEdge e : linkedList)
                {
                    // add the edge's cost to the tripCost
                    tripCost+=e.cost();
                    // print out both sides of the edge and its cost
                    int startCity = e.either();
                    System.out.println(cityNames[startCity] + " ✈︎✈︎✈︎✈︎ " + cityNames[e.other(startCity)] + " cost $" + e.cost());
                }
                // print out the tripCost
                System.out.println("Total trip cost: $" + tripCost + "\n");
            }
        }
    }

    private void addRoute()
    {
        boolean badInput = true; //for checking input validity
        int start = -1;
        int end = -1;
        double newCost = -1;
        int newdistance = -1;
        while(badInput){
            System.out.print("Please enter starting city (1-" + cityNames.length + "): ");
            start = Integer.parseInt(scan.nextLine())-1;
            System.out.print("Please enter ending city (1-" + cityNames.length + "): ");
            end = Integer.parseInt(scan.nextLine())-1;
            System.out.print("Please enter the cost between starting city and ending city: ");
            newCost = Double.parseDouble(scan.nextLine());
            System.out.print("Please enter the distance between starting city and ending city: ");
            newdistance = Integer.parseInt(scan.nextLine());
            //deal with the invalid input: number out of required bound, equal number, or a route is already existed between the two cities
            if((start>cityNames.length-1 || start<0) || (end>cityNames.length-1 || end<0) || start==end)
            {
                badInput = true;
                System.out.println("Invalid input. Enter again.");  
            }
            else if(routeExisted(start, end)) //check if a route is already existed between the two cities
            {
                badInput = true;
                System.out.println("The two cities you chose has a route already. Enter again.");
            }
            else
            {
                badInput = false;
            }
        }

        //if we go to here, we are ready to add the route
        G.addEdge(new DirectedEdge(start, end, newCost, newdistance));
        G.addEdge(new DirectedEdge(end, start, newCost, newdistance));

        System.out.println("New route added successfully!");
    }

    private boolean routeExisted(int start, int end)        //this method is to check if there is aleady a route between two specified cities
    {
        boolean exist = false;                              //flag

        outerloop:
        for(DirectedEdge e : G.adj(start))                  //iterate thru all edges of start city
        {
            if(e.to()==end)                                 //if there is an edge point to end city
            {
                for(DirectedEdge vEdge : G.adj(end))        //iterate thru all edges of end city 
                {
                    if((vEdge.to()==start))                 //if there is an edge point to start city
                    {
                        exist = true;                       //there is aleady a route between two specified cities
                        break outerloop;                    //stop to check 
                    }
                     
                }
            }
        }

        if(exist) return true;
        else return false;
    }

    private void removeRoute()                              //this has the similar logic as routeExisted() method
    {
        System.out.print("Please enter starting city (1-" + cityNames.length + "): ");
        int start = Integer.parseInt(scan.nextLine()) -1;
        System.out.print("Please enter ending city (1-" + cityNames.length + "): ");
        int end = Integer.parseInt(scan.nextLine()) -1;
        boolean removed = false;

        outerloop:
        for(DirectedEdge e : G.adj(start))
        {
            if(e.to()==end)
            {
                G.removeEdge(e);                            //when finding the edge, delete it
                for(DirectedEdge vEdge : G.adj(end))
                {
                    if((vEdge.to()==start))
                    {
                        G.removeEdge(vEdge);                //when finding the edge, delete it
                        removed = true;
                        G.decreaseNumEdges();
                        break outerloop;
                    }
                     
                }
            }
        }
        
        if(removed)
            System.out.println("Route has been successfully removed.");
        else
            System.out.println("Route cannot be removed.");

    }

    private void quitAndSave() throws IOException           
    {
        FileWriter writer = new FileWriter(fileName);
        boolean[][] finished = new boolean[G.adj.length][G.adj.length];
        writer.write(cityNames.length + "\n");
        for(String city : cityNames)
        {
            writer.write(city + "\n");
        }

        for(int i = 0; i<G.adj.length; i++)
        {
            for(DirectedEdge e : G.adj(i))
            {
                int start = e.from() +1;                //the purpose of increment and decrement this two vars is because
                int end = e.to() + 1;                   //the user views city starting from 1 but city start from 0 in code
                                                        //from() is to get the city from one side of edge
                                                        //to() is getting the other side city(edge)
                if(!finished[start-1][end-1])
                {
                    finished[start-1][end-1]=true;
                    finished[end-1][start-1]=true;
                    writer.write(start + " " + end + " " + e.distance() + " " + e.cost() + "\n");
                }
            }
        }
        writer.close();

        System.out.println("Your routes are all successfully saved into your file. Goodbye!");
        System.exit(0);
    }

    private void findShortestPathThruACity()                        //extra credit
    {
        System.out.print("Please enter starting city (1-" + cityNames.length + "): ");
        int source = Integer.parseInt(scan.nextLine())-1;
        System.out.print("Please enter ending city (1-" + cityNames.length + "): ");
        int destination = Integer.parseInt(scan.nextLine())-1;
        System.out.print("Please enter the third city between starting city and ending city: ");
        int middle = Integer.parseInt(scan.nextLine())-1;

        airline.thruACityHelper(source, middle); 
        airline.thruACityHelper(middle, destination);
        System.out.println("The shortest distance from " + cityNames[source] + " to " + cityNames[destination] + " via " + cityNames[middle] + " is " + shortestDistanceViaThirdCity + " miles");
    }

    private void thruACityHelper(int source, int destination)       //extra credit
    {                                                     //this helper method is pretty similar to getShortestPathByDistance() method
                                                         //the idea is to find the shortest path from A to C 
                                                         //and from C to B and add the two distance together
        G.digkstras(source, destination, false);          //the 3rd parameter is a boolean; false means by distance
                                                          //true means by cost
        if(!G.marked[destination])                        //marked array check is there is an edge between
        {                                                 //source and destination
            System.out.println("There is no route from " + cityNames[source] + " to " + cityNames[destination]);
        }
        else
        {
            Stack<Integer> path = new Stack<>();          //build a stack to store all the possible cities between s and destination
            for(int x = destination; x!=source; x=G.edgeTo[x]) //go to a previous edge on shortest path each time
            {
                path.push(x);
            }

            shortestDistanceViaThirdCity+=G.distTo[destination];  //increment the total shortest distance from A to C via B

            int preVertex = source;                       //pop those cities and calculate the distance and print out
            System.out.print(cityNames[source] + " ");
            while (!path.empty()) {
                int v = path.pop();
                System.out.print(G.distTo[v] - G.distTo[preVertex] + " " + cityNames[v] + " ");
                preVertex = v;
            }
            System.out.println();
        }
    }


   ////////////////////////other helper classes//////////////////////////

   /**
   * The <tt>Digraph</tt> class represents an directed graph of vertices named 0
   * through v-1. It supports the following operations: add an edge to the graph,
   * iterate over all of edges leaving a vertex.
   */
    private class Digraph                                   //this is from lab10
    {
        private final int v;        //number of vertices
        private int e;              //number of edges
        private LinkedList<DirectedEdge>[] adj;
        private boolean[] marked;   // marked[v] = is there an s-v path
        private int[] edgeTo;       // edgeTo[v] = previous edge on shortest s-v path
        private int[] distTo;       // distTo[v] = number of edges shortest s-v path

        //Create an empty digraph with v vertices.
        public Digraph(int v)
        {
            if(v<0) throw new RuntimeException("Number of vertices must be nonnegative");
            this.v=v;
            this.e=0;
            @SuppressWarnings("unchecked")
            LinkedList<DirectedEdge>[] temp = (LinkedList<DirectedEdge>[]) new LinkedList[v];
            adj = temp;
            for(int i = 0; i<v; i++)
            {
                adj[i] = new LinkedList<DirectedEdge>();
            }
        }

        public void decreaseNumEdges(){
            e--;
        }

        public void removeEdge( DirectedEdge edge){
            adj[edge.from()].remove(edge);
            
        }

        //Add the edge e to this digraph.
        public void addEdge(DirectedEdge edge)
        {
            int from = edge.from();
            adj[from].add(edge);
            e++;
        }

        public int v()
        {
            return this.v;
        }

        /**
         * Return the edges leaving vertex v as an Iterable. To iterate over the edges
         * leaving vertex v, use foreach notation:
         * <tt>for (DirectedEdge e : graph.adj(v))</tt>.
         */
        public Iterable<DirectedEdge> adj(int v)
        {
            return adj[v];
        }

        /**
         * Return all edges in this graph as an Iterable. To iterate over the edges, use
         * foreach notation: <tt>for (Edge e : graph.edges())</tt>.
         */
        public Iterable<DirectedEdge> edges() {
            ArrayList<DirectedEdge> list = new ArrayList<DirectedEdge>();
            for (int v = 0; v < v; v++) {
                int selfLoops = 0;
                for (DirectedEdge e : adj(v)) {
                    if (e.other(v) > v) {
                        list.add(e);
                    }
                    // only add one copy of each self loop
                    else if (e.other(v) == v) {
                        if (selfLoops % 2 == 0)
                            list.add(e);
                        selfLoops++;
                    }
                }
            }
            return list;
        }

        public void bfs(int source)
        {
            marked = new boolean[this.v];
            distTo = new int[this.v];
            edgeTo = new int[this.v];

            Queue<Integer> q = new LinkedList<Integer>();
            for(int i = 0; i<v; i++)
            {
                distTo[i] = INFINITY;
                marked[i] = false;
            }
            distTo[source] = 0;
            marked[source] = true;
            q.add(source);

            while(!q.isEmpty())
            {
                int v = q.remove();
                for(DirectedEdge w : adj(v))
                {
                    if(!marked[w.to()])
                    {
                        edgeTo[w.to()] = v;
                        distTo[w.to()] = distTo[v]+1;
                        marked[w.to()] = true;
                        q.add(w.to());
                    }
                }
            }
        }

        public void digkstras(int source, int destination, boolean costNotDistance)
        {
            marked = new boolean[this.v];
            distTo = new int[this.v];
            edgeTo = new int[this.v];

            for(int i = 0; i<v; i++)
            {
                distTo[i] = INFINITY;
                marked[i] = false;
            }
            distTo[source] = 0;
            marked[source] = true;
            int nMarked = 1;

            int current = source;
            while(nMarked<this.v)
            {
                for(DirectedEdge w : adj(current))
                {
                    if(costNotDistance)             //calculate depends on cost
                    {
                        if(distTo[current] + w.cost() < distTo[w.to()])
                        {
                            edgeTo[w.to()] = current;
                            distTo[w.to()] = distTo[current] + (int)w.cost;
                        }
                    }
                    else                            //calculate depends on distance
                    {
                        if(distTo[current] + w.distance() < distTo[w.to()])
                        {
                            edgeTo[w.to()] = current;
                            distTo[w.to()] = distTo[current] + w.distance;
                        }
                    }
                }

                // Find the vertex with minimim path distance
                // This can be done more effiently using a priority queue!
                int min = INFINITY;
                current = -1;

                for(int i = 0; i<distTo.length; i++)
                {
                    if(marked[i]) continue;
                    if(distTo[i] < min)
                    {
                        min = distTo[i];
                        current = i;
                    }
                }

                if(current>= 0)
                {
                    marked[current] = true;
                    nMarked++;
                }
                else break;
            }

        }
    }

   //The <tt>DirectedEdge</tt> class represents an edge in an directed graph.
    private class DirectedEdge                              //this is from lab10
    {
        private int v;
        private int w;
        private double cost;
        private int distance;

        public DirectedEdge(int v, int w, double cost, int distance)
        {
            this.v=v;
            this.w=w;
            this.setCost(cost);;
            this.setDistance(distance);
        }

        public int from()
        {
            return v;
        }

        public int to()
        {
            return w;
        }

        // Set the cost of this trip
        public void setCost(double cost) 
        {
            this.cost = cost;
        }

        public void setDistance(int distance) 
        {
            this.distance = distance;
        }

        public int distance()
        {
            return distance;
        }

        public double cost()
        {
            return cost;
        }

        public int either()
        {
            return v;
        }

        public int other(int vertex)
        {
            if (vertex == v)
                return w;
            else if (vertex == w)
                return v;
            else
                throw new RuntimeException("Illegal endpoint");
        }
    }

    public class PrimMSTTrace                               //this is from code handouts-Graphs
    {
        private DirectedEdge[] edgeTo;        // edgeTo[v] = shortest edge from tree vertex to non-tree vertex
        private double[] distTo;      // distTo[v] = weight of shortest such edge
        private boolean[] marked;     // marked[v] = true if v on tree, false otherwise
        private IndexMinPQ<Double> pq;

        public PrimMSTTrace(Digraph G) 
        {
            edgeTo = new DirectedEdge[G.v()];
            distTo = new double[G.v()];
            marked = new boolean[G.v()];
            pq = new IndexMinPQ<Double>(G.v());
            for (int v = 0; v < G.v(); v++) distTo[v] = Double.POSITIVE_INFINITY;

            for (int v = 0; v < G.v(); v++)      // run from each vertex to find
                if (!marked[v]) prim(G, v);      // minimum spanning forest

            // check optimality conditions
            assert check(G);
        }

        // run Prim's algorithm in graph G, starting from vertex s
        private void prim(Digraph G, int s) 
        {
            distTo[s] = 0.0;
            pq.insert(s, distTo[s]);
            //showPQ(pq);
            while (!pq.isEmpty()) 
            {
                int v = pq.delMin();
                //System.out.println("	Next Vertex (Weight): " + v + " (" + distTo[v] + ")");
                scan(G, v);
                //showPQ(pq);
            }
        }

        // scan vertex v
        private void scan(Digraph G, int v) 
        {
            marked[v] = true;
            //System.out.println("	Checking neighbors of " + v);
            for (DirectedEdge e : G.adj(v)) 
            {
                int w = e.other(v);
                //System.out.print("		Neighbor " + w);
                if (marked[w]) 
                {
                    //System.out.println(" is in the tree ");
                    continue;         // v-w is obsolete edge
                }
                if (e.distance() < distTo[w]) 
                {
                    //System.out.print(" OLD distance: " + distTo[w]);
                    distTo[w] = e.distance();
                    edgeTo[w] = e;
                    //System.out.println(" NEW distance: " + distTo[w]);
                    if (pq.contains(w)) 
                    {
                        pq.change(w, distTo[w]);
                    }
                    else              
                    {
                        pq.insert(w, distTo[w]);
                    }
                }
            }
        }
    
        // return iterator of edges in MST
        public Iterable<DirectedEdge> edges() {
            ArrayList<DirectedEdge> mst = new ArrayList<DirectedEdge>();
            for (int v = 0; v < edgeTo.length; v++) {
                DirectedEdge e = edgeTo[v];
                if (e != null) {
                    mst.add(e);
                }
            }
            return mst;
        }

        // return weight of MST
        public double weight() {
            double weight = 0.0;
            for (DirectedEdge e : edges())
                weight += e.distance();
            return weight;
        }

        // check optimality conditions (takes time proportional to E V lg* V)
        private boolean check(Digraph G) 
        {
            // check weight
            double weight = 0.0;
            for (DirectedEdge e : edges()) {
                weight += e.distance();
            }
            double EPSILON = 1E-12;
            if (Math.abs(weight - weight()) > EPSILON) {
                System.err.printf("Weight of edges does not equal weight(): %f vs. %f\n", weight, weight());
                return false;
            }

            // check that it is acyclic
            UF uf = new UF(G.v());
            for (DirectedEdge e : edges()) {
                int v = e.either(), w = e.other(v);
                if (uf.connected(v, w)) {
                    System.err.println("Not a forest");
                    return false;
                }
                uf.union(v, w);
            }

            // check that it is a spanning forest
            for (DirectedEdge e : edges()) {
                int v = e.either(), w = e.other(v);
                if (!uf.connected(v, w)) {
                    System.err.println("Not a spanning forest");
                    return false;
                }
            }

            // check that it is a minimal spanning forest (cut optimality conditions)
            for (DirectedEdge e : edges()) {
                int v = e.either(), w = e.other(v);

                // all edges in MST except e
                uf = new UF(G.v());
                for (DirectedEdge f : edges()) {
                    int x = f.either(), y = f.other(x);
                    if (f != e) uf.union(x, y);
                }

                // check that e is min weight edge in crossing cut
                for (DirectedEdge f : G.edges()) {
                    int x = f.either(), y = f.other(x);
                    if (!uf.connected(x, y)) {
                        if (f.distance() < e.distance()) {
                            System.err.println("Edge " + f + " violates cut optimality conditions");
                            return false;
                        }
                    }
                }

            }

            return true;
        }


}

    public class IndexMinPQ<Key extends Comparable<Key>> implements Iterable<Integer> ////this is from code handouts-Graphs
    {
        private int N;           // number of elements on PQ
        private int[] pq;        // binary heap using 1-based indexing
        private int[] qp;        // inverse of pq - qp[pq[i]] = pq[qp[i]] = i
        private Key[] keys;      // keys[i] = priority of i

        public IndexMinPQ(int NMAX) {
            keys = (Key[]) new Comparable[NMAX + 1];    // make this of length NMAX??
            pq   = new int[NMAX + 1];
            qp   = new int[NMAX + 1];                   // make this of length NMAX??
            for (int i = 0; i <= NMAX; i++) qp[i] = -1;
        }

        // is the priority queue empty?
        public boolean isEmpty() { return N == 0; }

        // is k an index on the priority queue?
        public boolean contains(int k) {
            return qp[k] != -1;
        }

        // number of keys in the priority queue
        public int size() {
            return N;
        }

        // associate key with index k
        public void insert(int k, Key key) {
            if (contains(k)) throw new RuntimeException("item is already in pq");
            N++;
            qp[k] = N;
            pq[N] = k;
            keys[k] = key;
            swim(N);
        }

        // return the index associated with a minimal key
        public int min() { 
            if (N == 0) throw new RuntimeException("Priority queue underflow");
            return pq[1];        
        }

        // return a minimal key
        public Key minKey() { 
            if (N == 0) throw new RuntimeException("Priority queue underflow");
            return keys[pq[1]];        
        }

        // delete a minimal key and returns its associated index
        public int delMin() { 
            if (N == 0) throw new RuntimeException("Priority queue underflow");
            int min = pq[1];        
            exch(1, N--); 
            sink(1);
            qp[min] = -1;            // delete
            keys[pq[N+1]] = null;    // to help with garbage collection
            pq[N+1] = -1;            // not needed
            return min; 
        }

    /*
        // change key associated with index k; insert if index k is not present
        public void put(int k, Key key) {
            if (!contains(k)) insert(k, key);
            else changeKey(k, key);
        }
        // return key associated with index k
        public Key get(int k) {
            if (!contains(k)) throw new RuntimeException("item is not in pq");
            else return keys[pq[k]];
        }
    */

        // change the key associated with index k
        public void change(int k, Key key) {
            if (!contains(k)) throw new RuntimeException("item is not in pq");
            keys[k] = key;
            swim(qp[k]);
            sink(qp[k]);
        }

        // decrease the key associated with index k
        public void decrease(int k, Key key) {
            if (!contains(k)) throw new RuntimeException("item is not in pq");
            if (keys[k].compareTo(key) <= 0) throw new RuntimeException("illegal decrease");
            keys[k] = key;
            swim(qp[k]);
        }

        // increase the key associated with index k
        public void increase(int k, Key key) {
            if (!contains(k)) throw new RuntimeException("item is not in pq");
            if (keys[k].compareTo(key) >= 0) throw new RuntimeException("illegal decrease");
            keys[k] = key;
            sink(qp[k]);
        }

    /**************************************************************
        * General helper functions
        **************************************************************/
        private boolean greater(int i, int j) {
            return keys[pq[i]].compareTo(keys[pq[j]]) > 0;
        }

        private void exch(int i, int j) {
            int swap = pq[i]; pq[i] = pq[j]; pq[j] = swap;
            qp[pq[i]] = i; qp[pq[j]] = j;
        }

    /**************************************************************
        * Heap helper functions
        **************************************************************/
        private void swim(int k)  {
            while (k > 1 && greater(k/2, k)) {
                exch(k, k/2);
                k = k/2;
            }
        }

        private void sink(int k) {
            while (2*k <= N) {
                int j = 2*k;
                if (j < N && greater(j, j+1)) j++;
                if (!greater(k, j)) break;
                exch(k, j);
                k = j;
            }
        }


    /***********************************************************************
        * Iterators
        **********************************************************************/

    /**
         * Return an iterator that iterates over all of the elements on the
         * priority queue in ascending order.
         * <p>
         * The iterator doesn't implement <tt>remove()</tt> since it's optional.
         */
        public Iterator<Integer> iterator() { return new HeapIterator(); }

        private class HeapIterator implements Iterator<Integer> {
            // create a new pq
            private IndexMinPQ<Key> copy;

            // add all elements to copy of heap
            // takes linear time since already in heap order so no keys move
            public HeapIterator() {
                copy = new IndexMinPQ<Key>(pq.length - 1);
                for (int i = 1; i <= N; i++)
                    copy.insert(pq[i], keys[pq[i]]);
            }

            public boolean hasNext()  { return !copy.isEmpty();                     }
            public void remove()      { throw new UnsupportedOperationException();  }

            public Integer next() {
                if (!hasNext()) throw new NoSuchElementException();
                return copy.delMin();
            }
        }
    }

    /**
 *  The <tt>UF</tt> class represents a union-find data data structure.
 *  It supports the <em>union</em> and <em>find</em>
 *  operations, along with a method for determining the number of
 *  disjoint sets.
 *  <p>
 *  This implementation uses weighted quick union.
 *  Creating a data structure with N objects takes linear time.
 *  Afterwards, all operations are logarithmic worst-case time.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/15uf">Section 1.5</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 */
    public class UF {                                       //this is from code handouts-Graphs
        private int[] id;    // id[i] = parent of i
        private int[] sz;    // sz[i] = number of objects in subtree rooted at i
        private int count;   // number of components

    /**
         * Create an empty union find data structure with N isolated sets.
         */
        public UF(int v) {
            count = v;
            id = new int[v];
            sz = new int[v];
            for (int i = 0; i < v; i++) {
                id[i] = i;
                sz[i] = 1;
            }
        }

    /**
         * Return the id of component corresponding to object p.
         */
        public int find(int p) {
            int length = id.length;
            if (p < 0 || p >= length) {
                throw new IndexOutOfBoundsException("Cannot access " + p);
            }

            while (p != id[p])
                p = id[p];
            return p;
        }

    /**
         * Return the number of disjoint sets.
         */
        public int count() {
            return count;
        }

        public int[] parentArr()
        {
            return this.id;
        }

        public int[] sizeArr()
        {
            return this.sz;
        }
    
    /**
         * Are objects p and q in the same set?
         */
        public boolean connected(int p, int q) {
            return find(p) == find(q);
        }

    
    /**
         * Replace sets containing p and q with their union.
         */
        public void union(int p, int q) {
            int i = find(p);
            int j = find(q);
            if (i == j) return;

            // make smaller root point to larger one
            if   (sz[i] < sz[j]) { id[i] = j; sz[j] += sz[i]; }
            else                 { id[j] = i; sz[i] += sz[j]; }
            count--;
        }

    }

    public class DepthFirstAmountPaths                      //this is from code handouts-Graphs
    {
        private final double amount; // max amount
        ArrayList<LinkedList<DirectedEdge>> results = new ArrayList<LinkedList<DirectedEdge>>();
    
        public DepthFirstAmountPaths(Digraph G, int s, double amount) 
        {
            this.amount = amount;
            dfs(G, s, new LinkedList<DirectedEdge>(), 0, new boolean[G.v()]);
        }

        // depth first search from v
        private void dfs(Digraph G, int v, LinkedList<DirectedEdge> currentPath, double currAmount, boolean[] marked) 
        {
            marked[v] = true;
    
            int w;
            for (DirectedEdge e : G.adj(v)) {
                w = e.to();
    
                // If the edge's destination has not yet been visited
                if (!marked[w]) {
                    double newAmount = currAmount + e.cost();
    
                    // If the cost of the trip does not cause us to bust
                    if (newAmount <= amount) {
                        LinkedList<DirectedEdge> newPath = cloneList(currentPath);
                        newPath.add(e);
                        results.add(newPath);
                        boolean[] newMarked = cloneArr(marked);
    
                        // Recurse to destination with stack frame-specific variables
                        dfs(G, w, newPath, newAmount, newMarked);
                    }
                }
            }
        }
    
        public ArrayList<LinkedList<DirectedEdge>> getResults() {
            return results;
        }
    
        private LinkedList<DirectedEdge> cloneList(LinkedList<DirectedEdge> oldList) {
            LinkedList<DirectedEdge> newList = new LinkedList<DirectedEdge>();
    
            for (DirectedEdge e : oldList) {
                newList.add(e);
            }
    
            return newList;
        }
    
        private boolean[] cloneArr(boolean[] oldArr) {
            boolean[] newArr = new boolean[oldArr.length];
    
            for (int i = 0; i < oldArr.length; i++) {
                newArr[i] = oldArr[i];
            }
    
            return newArr;
        }
    
    }
}
