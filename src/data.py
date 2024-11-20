import numpy as np


class Data:
    def __init__(self, filename: str):
        self.factor = 0.39365
        self.way_cost = 1000 * self.factor
        self.edge_len = np.array([8.500, 20.870])
        self.primary_cost = 27.7 * self.factor
        self.secondary_cost = 10 * self.factor
        self.capacity = [6, 8]
        self.price = {6: 240 * self.factor, 8: 280 * self.factor}

        self.location = []
        self.belong = []
        self.district_size = 109

        with open(filename, 'r') as f:
            self._initialize_from_file(f)

        self._calculate_way_attributes()
        self._process_graph()
        self._calculate_distances()
        self._initialize_way_set()

    def _initialize_from_file(self, file):
        line = file.readline().split()
        self.row_num = int(line[0])
        self.col_num = int(line[1])
        self.node_num = int(line[2])
        self.district_num = self.node_num // self.district_size
        self.district = [[] for _ in range(self.district_num)]
        self.graph = np.zeros((self.row_num, self.col_num), dtype=int)

        for cnt, line in enumerate(file.readlines()):
            x, y, r = map(int, line.split())
            self.graph[x, y] = r
            self.belong.append(r - 1)
            self.district[r - 1].append(cnt)
            self.location.append(np.array([x, y]))

    def _calculate_way_attributes(self):
        self.way_num = self.col_num - 1
        self.way_len = [
            max(np.count_nonzero(self.graph[:, i]), np.count_nonzero(self.graph[:, i + 1]))
            for i in range(self.way_num)
        ]
        self.way_cover = [
            [np.any(self.graph[:, i:i + 2] - 1 == k) for k in range(self.district_num)]
            for i in range(self.way_num)
        ]

    def _process_graph(self):
        self.graph[self.graph == 0] = -1
        for i, (x, y) in enumerate(self.location):
            self.graph[x, y] = i

    def _calculate_distances(self):
        self.distance = np.array([
            [
                sum(abs(self.location[i] - self.location[j]) * self.edge_len)
                for j in range(self.node_num)
            ]
            for i in range(self.node_num)
        ])

    def _initialize_way_set(self):
        self.way_set = [
            np.concatenate((self.graph[:, i], self.graph[:, i + 1]))
            for i in range(self.way_num)
        ]
        self.way_set = [x[x >= 0] for x in self.way_set]
