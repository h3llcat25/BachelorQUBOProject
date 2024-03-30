class TieQuboStruct:
    def __init__(self, bin_var_dict=None, objective_propagation=None, objective_damage=None, graph=None, target_nodes=None):
        self.bin_var_dict = bin_var_dict if bin_var_dict is not None else {}
        self.objective_propagation = objective_propagation if objective_propagation is not None else []
        self.objective_damage = objective_damage if objective_damage is not None else []
        self.graph = graph if graph is not None else []
        self.target_nodes = target_nodes if target_nodes is not None else []

    def set_bin_var_dict(self, bin_var_dict):
        self.bin_var_dict = bin_var_dict

    def get_dict_objectives_graph_targets(self):
        return self.bin_var_dict, self.objective_propagation, self.objective_damage, self.graph, self.target_nodes

    def get_dict_objectives_graph(self):
        return self.bin_var_dict, self.objective_propagation, self.objective_damage, self.graph

    def has_seed(self):
        return self.seed is not None
