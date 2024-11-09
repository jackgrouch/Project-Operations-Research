# Result Description

Each TXT file named "district_<span>$a$</span>-<span>$b$</span>.txt" contains the optimal solution for the district_<span>$a$</span>-<span>$b$</span> instance, solved using our exact solution approach. The file structure is as follows:

* The first line reports the optimal value in USD.
* The second line, labeled "Service ways", lists the indices of service ways to be built.
* The third line, labeled "Inverters", lists the indices of PVAs designated for inverter installations.
* After the fourth line, labeled "Cable routing", each subsequent line reports the cable routing solution for each district based on the following `dictionary` structure:

  * `key`: the index of PVA where the combiner box is located.
  * `value`: the indices of PVAs connected to the combiner box.
