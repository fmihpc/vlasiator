class MaskEveryNthCell(SimplePythonExpression):
    def __init__(self):
        SimplePythonExpression.__init__(self)
        self.name = "MaskEveryNthCell"
        self.description = "Set a cell value to 0 if it should be masked, 1 if it should not."
        self.output_is_point_var  = False
        self.output_dimension = 1
        # parameters to change
        xCells = 666
        yCells = 583
        xStep = 50
        yStep = 50
        xList = range(xStep, xCells+1, xStep)
        yList = range(yStep, yCells+1, yStep)
        self.cellList = []
        for x in xList:
            for y in yList:
                self.cellList.append((y-1)*xCells + x)

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
 
py_filter = MaskEveryNthCell
