object elim {
	val node1 = new Node(1)                   //> node1  : Node = Node@71971d92
	val node2 = new Node(2)                   //> node2  : Node = Node@1a21a14e
	val node3 = new Node(3)                   //> node3  : Node = Node@15e51549
	val node4 = new Node(4)                   //> node4  : Node = Node@72d876d9
	val edge1 = ((node1, node2))              //> edge1  : (Node, Node) = (Node@71971d92,Node@1a21a14e)
	val edge2 = ((node2, node3))              //> edge2  : (Node, Node) = (Node@1a21a14e,Node@15e51549)
	val edge3 = ((node3, node4))              //> edge3  : (Node, Node) = (Node@15e51549,Node@72d876d9)
	val G = new Graph(List(node1, node2, node3, node4), List(edge1, edge2, edge3))
                                                  //> G  : Graph = Graph@64bf091d
	val I = List(4,3,2,1) // eliminate ordering
                                                  //> I  : List[Int] = List(4, 3, 2, 1)
	
	def eliminate(G: Graph, E: List[Int], F: List[Int]) = {
		val (I, probList) = initialize(G, F)
		val evidenceList = evidence(E)
		val prob = update(G, I, probList:::evidenceList).head.values
		normalize(F, prob)
	}                                         //> eliminate: (G: Graph, E: List[Int], F: List[Int])Array[Double]
	
	def initialize(G: Graph, F: List[Int]) = {
		val I = List(4,3,2,1)
		val activeList:List[Potential] = G.nodes map ( _.localProb )
		(I, activeList)
	}                                         //> initialize: (G: Graph, F: List[Int])(List[Int], List[Potential])
	
	def evidence(E: List[Int]) = {
		E map { Evidence(_, 1) }
	}                                         //> evidence: (E: List[Int])List[Evidence]
	
	def update(G: Graph, I: List[Int], activeList: List[Potential]) = {
		def updateIter(I: List[Int], activeList: List[Potential]): List[Potential] = {
			println("Active List => " + activeList)
			
			def find(list: List[Potential], i: Int): (List[Potential], List[Potential], Int, Int) = {
				list match {
					case List() => (List(), List(), 0, -1)
					case head::tail => {
						val (phi, newList, si, evidence) = find(tail, i)
						head match {
							case Prob(x) =>
								if (x == i)	(head::phi, newList, si, evidence)
								else (phi, head::newList, si, evidence)
							case CondProb(x, y) =>
								if (x == i) (head::phi, newList, y, evidence)
								else (phi, head::newList, si, evidence)
							case Evidence(x, y) =>
								if (x == i) (phi, newList, si, y)
								else (phi, head::newList, si, evidence)
							case M(x, y) =>
								if (x == i) (head::phi, newList, si, evidence)
								else (phi, head::newList, si, evidence)
							case _ => (phi, head::newList, si, evidence)
						}
					}
				}
			}
			
			def sum (p1: Array[Double], p2: Array[Double]) = {
				if (p2.size == 4) {
					val t00 = p1(0) * p2(0)
					val t01 = p1(0) * p2(1)
					val t10 = p1(1) * p2(2)
					val t11 = p1(1) * p2(3)
					Array(t00 + t10, t01 + t11)
				} else {
					val t0 = p1(0) * p2(0)
					val t1 = p1(1) * p2(1)
					Array(t0, t1)
				}
			}
			
			I match {
				case List() => activeList
				case i::tail => {
					val (phi, newList, si, evidence) = find (activeList, i)
					val prods = phi map ( _.values )
					if (evidence > 0) {
						val t = prods.head
						val m = Array(t(2*evidence), t(2*evidence + 1))
						println("m%d(x%d) => ".format(i, si) + m.toList)
						updateIter(tail, M(si, m)::newList)
					}
					else {
						val sortedProds = prods.sortWith (_.size < _.size)
						val m = sum(sortedProds.head, sortedProds.last)
						println("m%d(x%d) => ".format(i, si) + m.toList)
						updateIter(tail, M(si, m)::newList)
					}
				}
			}
		}
		
		updateIter(I, activeList)
	}                                         //> update: (G: Graph, I: List[Int], activeList: List[Potential])List[Potential
                                                  //| ]
	
	def normalize(F: List[Int], prob: Array[Double]) = {
		val sum = prob(0) + prob(1)
		Array(prob(0) / sum, prob(1) / sum)
	}                                         //> normalize: (F: List[Int], prob: Array[Double])Array[Double]
	
	val px1x4 = eliminate(G, List(4), List(1))//> Active List => List(Prob(1), CondProb(2,1), CondProb(3,2), CondProb(4,3), E
                                                  //| vidence(4,1))
                                                  //| m4(x3) => List(0.4, 0.8)
                                                  //| Active List => List(M(3,[D@29c82eaa), Prob(1), CondProb(2,1), CondProb(3,2)
                                                  //| )
                                                  //| m3(x2) => List(0.56, 0.7200000000000002)
                                                  //| Active List => List(M(2,[D@5cc75d38), Prob(1), CondProb(2,1))
                                                  //| m2(x1) => List(0.6240000000000001, 0.6880000000000002)
                                                  //| Active List => List(M(1,[D@24b9371e), Prob(1))
                                                  //| m1(x0) => List(0.31200000000000006, 0.3440000000000001)
                                                  //| Active List => List(M(0,[D@7bd1a567))
                                                  //| px1x4  : Array[Double] = Array(0.47560975609756095, 0.524390243902439)
	px1x4(0)                                  //> res0: Double = 0.47560975609756095
}

abstract class Potential{
	val values: Array[Double]
}
case class Prob (x: Int) extends Potential {
	val values = Array(0.5, 0.5)
}
case class CondProb (x: Int, y: Int) extends Potential {
  val values = Array(0.6, 0.2, 0.4, 0.8)
}
case class Evidence (x: Int, y: Int) extends Potential {
	val values = Array(0.0)
}
case class M (x: Int, val values: Array[Double]) extends Potential


class Node (val idx: Int) {
	val localProb = if (idx == 1) Prob(idx) else CondProb(idx, idx - 1)
}

class Graph (val nodes: List[Node], val edges: List[(Node, Node)]) {
	def node(idx: Int) = {
		var n: Node = null
		for (node <- nodes) {
		  if (node.idx == idx)
		  	n = node
		}
		n
	}
}