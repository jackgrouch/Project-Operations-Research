# Instance Description

Each PNG file "district_<span>$a$</span>-<span>$b$</span>.png" in [instances_picture](instances_picture) directory illustrates the placement of PVAs, where the partitioning result is obtained by our Heuristic Search algorithm, and the blue dotted line represents the induced inner cut edges.

Each TXT file "district_<span>$a$</span>-<span>$b$</span>.txt" in [instances](instances) directory provides information about the location of PVAs, with each line reporting the PVA's index, x-coordinate, y-coordinate, and the district to which the PVA is assigned.

The length parameters are provided in meters. The distance between two PVAs is measured by the Manhattan distance (L1 norm).  The potential service ways can be built by vertically traversing the entire farm between any two consecutive PVA columns. The length of a service way is calculated as $8.500 \times max(n^{1},n^{2})$, where $n^{1}$ and $n^{2}$ are respectively the number of PVAs located along the two adjacent columns of this service way.

The cost parameters are provided in USD. The cabling cost for connecting a combiner box to the inverter is $c_1 = 10.904$ per meter and that for connecting a PVA to a combiner box is $c_2 = 3.937$ per meter. The cost for building a service way is $46.312$
per meter. Two types of combiner box ($T=2$), each with capacity $p_1=6$ and $p_2 =8$, are used in these plants and their purchasing costs $g_1 = 94.476$ and $g_2 = 110.222$, respectively.
