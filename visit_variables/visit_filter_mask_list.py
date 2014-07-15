class MaskList(SimplePythonExpression):
    def __init__(self):
        SimplePythonExpression.__init__(self)
        self.name = "MaskList"
        self.description = "Set a cell value to 0 if it should be masked, 1 if it should not."
        self.output_is_point_var  = False
        self.output_dimension = 1
        # parameters to change
        self.cellList = [71680,85390,103324,148290,148436,148456,148462,157260,175360]

    def derive_variable(self,ds_in,domain_id):
        # ds_in is a vtk dataset, we want
        # to create and return a new vtkDataArray
        # that contains the mask pattern
        cell_ids = ds_in.GetCellData().GetArray(self.input_var_names[0])
        res = vtk.vtkFloatArray()
        ncells = ds_in.GetNumberOfCells()
        res.SetNumberOfComponents(1)
        res.SetNumberOfTuples(ncells)
        for i in xrange(ncells):
            cell = ds_in.GetCell(i)
            cellid = cell_ids.GetTuple1(i)
            val = 1.0
            if cellid in self.cellList:
                val = 0.0
            res.SetTuple1(i,val)
        return res
 
py_filter = MaskList
