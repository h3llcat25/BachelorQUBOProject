class TieGraphProblem:
    def __init__(self, bin_var_dict=None, list1=None, list2=None, list3=None, list4=None, seed=None):
        self.bin_var_dict = bin_var_dict if bin_var_dict is not None else {}
        self.marked_c_nodes = list1 if list1 is not None else []
        self.c_nodes = list2 if list2 is not None else []
        self.e_nodes = list3 if list3 is not None else []
        self.r_nodes = list4 if list4 is not None else []
        self.seed = seed

        self.abstrct_dict = None
        self.equation = ''
        self.react_pena_term = ''
        self.extra_var_index_r = 1
        self.abstract_term = ''

    def set_bin_var_dict(self, bin_var_dict):
        self.bin_var_dict = bin_var_dict

    def get_bin_var_dict(self):
        return self.bin_var_dict

    def set_marked_c_nodes(self, marked_c_nodes):
        self.marked_c_nodes = marked_c_nodes

    def get_marked_c_nodes(self):
        return self.marked_c_nodes

    def set_c_nodes(self, c_nodes):
        self.c_nodes = c_nodes

    def set_r_nodes(self, r_nodes):
        self.r_nodes = r_nodes

    def set_e_nodes(self, e_nodes):
        self.e_nodes = e_nodes

    def set_seed(self, seed):
        self.seed = seed

    def get_seed(self):
        return self.seed

    def set_equation(self, equation):
        self.equation = equation

    def get_equation(self):
        return self.equation

    def set_react_pena_term(self, react_pena_term):
        self.react_pena_term = react_pena_term

    def get_react_pena_term(self):
        return self.react_pena_term

    def set_extra_var_index_r(self, extra_var_index_r):
        self.extra_var_index_r = extra_var_index_r

    def get_extra_var_index_r(self):
        return self.extra_var_index_r

    def get_dict_and_lists(self):
        return self.bin_var_dict, self.marked_c_nodes, self.c_nodes, self.r_nodes, self.e_nodes

    def has_seed(self):
        return self.seed is not None

    def set_abstrct_dict(self, abstrct_dict):
        self.abstrct_dict = abstrct_dict

    def get_abstrct_dict(self):
        return self.abstrct_dict

    def set_abstract_term(self, abstract_term):
        self.abstract_term = abstract_term

    def get_abstract_term(self):
        return self.abstract_term