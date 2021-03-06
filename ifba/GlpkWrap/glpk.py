#!/usr/bin/env python
# encoding: utf-8
"""
glpk.py

Created by Nikolaus Sonnenschein on 2008-02-12.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import os
import sys
from ifba.glpki.glpki import *
from ifba.GlpkWrap import util, exceptions
from ifba.general.util import randomString

class glpk(object):
    """docstring for GlpkWrapper"""
    def __init__(self, lp):
        """Initialization of class attributes.
        smcp: simplex parameters. The default message level is that all output
        is allowed. The lp gets also indexed so that string ids of rows and
        columns can be used."""
        glp_term_out(GLP_OFF)
        self.lp = lp
        self.smcp = glp_smcp()
        self.iocp = glp_iocp()
        glp_init_smcp(self.smcp)
        glp_init_iocp(self.iocp)
        self.smcp.msg_lev = GLP_MSG_OFF
        self.iocp.msg_lev = GLP_MSG_OFF
        self.smcp.presolve = GLP_OFF
        self.iocp.presolve = GLP_ON
        glp_create_index(self.lp)
        glp_scale_prob(self.lp, GLP_SF_AUTO)
        self.history = []
        self.errDict = {
        GLP_UNDEF:"/* solution is undefined */",
        GLP_UNDEF:"/* solution is undefined */",
        GLP_FEAS:"/* solution is feasible */",
        GLP_INFEAS:"/* solution is infeasible */",
        GLP_NOFEAS:"/* no feasible solution exists */",
        GLP_OPT:"/* solution is optimal */",
        GLP_UNBND:"/* solution is unbounded */"
        }
    
    def __str__(self):
        """Prints information about the glpk object and it's current state."""
        optMappint = {1:'Minimize', 2:'Maximize'}
        verbMapping = {0:'All Messages Off', 1:'Error and Warnings only',
        2:'Normal Output', 3:'All Messages On'}
        presolveMapping = {0:'Off', 1:'On'}
        info = (
        """Some info:
        The current glpk version: %s
        The current Optimiztion direction: %s
        The current objective value: %f
        The number of Columns: %d
        The number of Rows: %d
        The number of Integer Columns: %d
        The number of Binary Columns: %d
        The number of undoable steps: %d
        The current verbosity level: %s
        Presolver: %s
        """ % (glp_version() ,optMappint[self.getOptFlag()], self.getObjVal(),
        self.getNumCols(), self.getNumRows(), glp_get_num_int(self.lp), glp_get_num_bin(self.lp),
        len(self.history), verbMapping[self.smcp.msg_lev],
        presolveMapping[self.smcp.presolve])
        )
        return str(info)
    
    def cplex(self):
        """docstring for __repr__"""
        rndFileStr = randomString(10)
        util.WriteCplex(self, rndFileStr)
        f = open(rndFileStr).read()
        os.remove(rndFileStr)
        return f
    
    def writeSolution(self, file="/dev/stdout"):
        glp_write_sol(self.lp, file)
        
    def undo(self):
        """Reverts the last modification."""
        try:
            self.history[0].undo()
            self.history.pop(0)
        except: IndexError, "no undo available"
    
    def initialize(self):
        """Reverts all modifications."""
        for undoObj in self.history:
            undoObj.undo()
        self.history = []
        
    def eraseHistory(self):
        """Erases all stored modifications."""
        self.history = []
    
    def __copy__(self):
        """Special method defined for usage with copy module.
        copy.copy(lp) yields a copy of lp.
        """
        newProb = glp_create_prob()
        glp_copy_prob(newProb, self.lp, GLP_ON)
        return glpk(newProb)
    
    def __getstate__(self):
        """How to pickle.""" #TODO implement without temporary file workaround. Must use some pipes and stdout
        return self.__repr__()
    
    def __setstate__(self, stuff):
        """How to unpickle."""
        rndFileStr = randomString(10)
        open(rndFileStr, 'w').write(stuff)
        self.__init__(util.ImportCplex(rndFileStr))
        os.remove(rndFileStr)
    
    def toggleVerbosity(self):
        """Toggles the verbosity level of glpk.
        The toggle switches between all messages allowed to no messages
        allowed."""
        # return glp_term_out(GLP_ON)
        if self.smcp.msg_lev == GLP_MSG_ON:
            self.smcp.msg_lev = GLP_MSG_OFF
        elif self.smcp.msg_lev == GLP_MSG_OFF:
            self.smcp.msg_lev = GLP_MSG_ON
        if self.iocp.msg_lev == GLP_MSG_ON:
            self.iocp.msg_lev = GLP_MSG_OFF
        elif self.iocp.msg_lev == GLP_MSG_OFF:
            self.iocp.msg_lev = GLP_MSG_ON
    
    def getNumCols(self):
        """Returns the number of columns in lp."""
        return glp_get_num_cols(self.lp)
    
    def getNumRows(self):
        """Returns the number of rows int lp."""
        return glp_get_num_rows(self.lp)
    
    def simplex(self):
        """Solves the lp using the glpk simplex function. Returns the return
        value of the glpk method."""
        self.smcp.meth = GLP_PRIMAL
        glp_simplex(self.lp, self.smcp)
        status = glp_get_status(self.lp)
        if status != GLP_OPT:
            if status == GLP_UNBND:
                raise exceptions.SolutionUnbounded, str(self.errDict[status])
            elif status == GLP_INFEAS:
                raise exceptions.SolutionInfeasible, str(self.errDict[status])
            elif status == GLP_NOFEAS:
                raise exceptions.NoFeasibleSolution, str(self.errDict[status])
            elif status == GLP_UNDEF:
                raise exceptions.SolutionUndefined, str(self.errDict[status])
            else:
                raise Exception, str(self.errDict[status])
        return status
    
    def exact(self):
        """Solves the lp with exact arithmetic using lpx_exact. Return the
        return value of the glpk function."""
        glp_exact(self.lp, self.smcp)
        status = glp_get_status(self.lp)
        if status != GLP_OPT:
            raise Exception, str(self.errDict[status])
        return status
    
    def interiorPoint(self):
        """Solves the lp using the glpk primal-dual interior-point method."""
        # FIXME Problems with numerical instabilities -> To large bounds could be the problem
        return glp_interior(self.lp, None)
        
    def intopt(self):
        """Solves a milp using glp_intopt."""
        glp_intopt(self.lp, self.iocp)
        status = glp_mip_status(self.lp)
        if status != GLP_OPT:
            raise Exception, str(self.errDict[status])
        return status
        
    def getObjVal(self):
        """Returns the current objective value."""
        return glp_get_obj_val(self.lp)

    def getInteriorObjVal(self):
        """Returns the current objective value."""
        return glp_get_obj_val(self.lp)
    
    def _checkColumnIndexValidity(self, index):
        numCols = self.getNumCols()
        if index < 1 or index > numCols:
            errorString = "Column index "+str(index)+" out of range.\n"
            errorString2 = "The valid range is: 1 - " + repr(numCols)
            raise IndexError, errorString + errorString2
    
    def _checkRowIndexValidity(self, index):
        numRows = self.getNumRows()
        if index < 1 or index > numRows:
            errorString = "Row index "+str(index)+" out of range.\n"
            errorString2 = "The valid range is: 1 - " + repr(numRows)
            raise IndexError, errorString + errorString2
    
    def _setColumnBound(self, index, lb, ub):
        """Sets a column bound.
        Determines the appropriate column type by itself."""
        self._checkColumnIndexValidity(index)
        if lb == '-inf' and ub == 'inf':
            glp_set_col_bnds(self.lp, index, GLP_FR, 0., 0.) # 0.'s are ignored
        elif lb == '-inf':
            glp_set_col_bnds(self.lp, index, GLP_UP, 0., ub) # 0. is ignored
        elif ub == 'inf':
            glp_set_col_bnds(self.lp, index, GLP_LO, lb, 0.) # 0. is ignored
        elif lb == ub:
            glp_set_col_bnds(self.lp, index, GLP_FX, lb, ub)
        elif lb != ub:
            glp_set_col_bnds(self.lp, index, GLP_DB, lb, ub)
        else:
            raise Exception, "Something is wrong with the provided bounds " + str(lb) + " " + str(ub)
    
    def _setRowBound(self, index, lb, ub):
        """Sets a column bound.
        Determines the appropriate column type by itself."""
        self._checkRowIndexValidity(index)
        if lb == '-inf' and ub == 'inf':
            glp_set_row_bnds(self.lp, index, GLP_FR, 0., 0.) # 0.'s are ignored
        elif lb == '-inf':
            glp_set_row_bnds(self.lp, index, GLP_UP, 0., ub) # 0. is ignored
        elif ub == 'inf':
            glp_set_row_bnds(self.lp, index, GLP_LO, lb, 0.) # 0. is ignored
        elif lb == ub:
            glp_set_row_bnds(self.lp, index, GLP_FX, lb, ub)
        elif lb != ub:
            glp_set_row_bnds(self.lp, index, GLP_DB, lb, ub)
        else:
            raise Exception, "Something is wrong with the provided bounds " + str(lb) + " " + str(ub)
    
    def _modifyColumnBoundsNoMemory(self, boundDict):
        """Modifies column bounds in lp without using the memory of lp."""
        for elem in boundDict:
            bnds = boundDict[elem]
            self._setColumnBound(elem, bnds[0], bnds[1])
    
    def modifyColumnBounds(self, boundDict):
        """Modifies column bounds of lp according to dictionary of bounds, e.g. {1:(0, 20), 44:(0, 20)}, to lp.
        Uses the memory of lp.
        """
        command = SetColumnBoundsCommand(self, boundDict)
        command.execute()
        
    def modifyRowBounds(self, boundDict):
        """Modifies row bounds of lp according to dictionary of bounds, e.g. {1:(0, 20), 44:(0, 20)}, to lp.
        Uses the memory of lp.
        """
        command = SetRowBoundsCommand(self, boundDict)
        command.execute()

    def _modifyRowBoundsNoMemory(self, boundDict):
        """Modifies row bounds in lp without using the memory of lp."""
        for elem in boundDict:
            bnds = boundDict[elem]
            self._setRowBound(elem, bnds[0], bnds[1])
    
    def getColumnBounds(self):
        """Return a dict of bounds."""
        boundsDict = dict()
        for elem in xrange(1, self.getNumCols() + 1):
            boundsDict[elem] = (glp_get_col_lb(self.lp, elem),
                glp_get_col_ub(self.lp, elem))
        return boundsDict

    def getRowBounds(self):
        """Return a dict of bounds."""
        boundsDict = dict()
        for elem in xrange(1, self.getNumRows() + 1):
            boundsDict[elem] = (glp_get_row_lb(self.lp, elem),
                glp_get_row_ub(self.lp, elem))
        return boundsDict
    
    def getColumnIDs(self):
        """Return a tuple of columnIDs"""
        num = self.getNumCols()
        return tuple([glp_get_col_name(self.lp, i) for i in range(1, num+1)])

    def getRowIDs(self):
        """Return a tuple of columnIDs"""
        num = self.getNumRows()
        return tuple([glp_get_row_name(self.lp, i) for i in range(1, num+1)])
    
    def translateColumnIndices(self, indices):
        """Translates a list of column numbers to the corresponding column
        ids."""
        for index in indices:
            self._checkColumnIndexValidity(index)
        return [glp_get_col_name(self.lp, i) for i in indices]
    
    def translateRowIndices(self, indices):
        """Translates a list of row numbers to the corresponding column ids."""
        for index in indices:
            self._checkRowIndexValidity(index)
        return [glp_get_row_name(self.lp, i) for i in indices]
    
    def translateColumnNames(self, names):
        """Translates a list of column ids to the corresponding column ordinal
        numbers."""
        ret = list()
        for name in names:
            index = glp_find_col(self.lp, name)
            if index != 0:
                ret.append(index)
            else:
                raise Exception, "Column id " + name + " cannot be found."
        return ret
    
    def translateRowNames(self, names):
        """Translates a list of row ids to the corresponding column ordinal
        numbers."""
        ret = list()
        for name in names:
            index = glp_find_row(self.lp, name)
            if index != 0:
                ret.append(index)
            else:
                raise Exception, "Row id " + name + " cannot be found."
        return ret
    
    def deleteColumns(self, colIndexes):
        """Deletes the columns specified in colIndexes."""
        for index in colIndexes:
            self._checkColumnIndexValidity(index)
        command = DeleteColumnsCommand(self, colIndexes)
        command.execute()
    
    def deleteRows(self, rowIndexes):
        """Deletes the rows specified in rowIndexes."""
        for index in rowIndexes:
            self._checkRowIndexValidity(index)
        command = DeleteRowsCommand(self, rowIndexes)
        command.execute()
    
    def addColumns(self, columns):
        """Adds new columns to the lp.
        columns is a dictionary of the form {columnID : (lb, ub, coeffList)}
        """
        command = AddColumnsCommand(self, columns)
        command.execute()
    
    def setColumnKinds(self, kindsDict):
        """Sets kinds of columns specified in kindsDict
        
        kindsDict has the fomr
        """
        command = SetColumnKindsCommand(self, kindsDict)
        command.execute()
    
    def getColumnKinds(self, colIndexes=None):
        columnKinds = dict()
        if not colIndexes:
            colIndexes = xrange(1, self.getNumCols() + 1)
        for index in colIndexes:
            self._checkColumnIndexValidity(index)
        for i in colIndexes:
            columnKinds[i] = glp_get_col_kind(self.lp, i)
        return columnKinds
    
    def addRows(self, rows):
        """Adds new rows to the lp.
        rows is a dictionary of the form {rowID : (lb, ub, coeffList)}
        """
        command = AddRowsCommand(self, rows)
        command.execute()
    
    def getColumnCoef(self, colIndex):
        """Returns a list of the column coefficients specified by colIndex."""
        self._checkColumnIndexValidity(colIndex)
        ia = intArray(self.getNumRows() + 1)
        da = doubleArray(self.getNumRows() + 1)
        numNonZero = glp_get_mat_col(self.lp, colIndex, ia, da)
        colCoefList = [0. for i in range(self.getNumRows())]
        for i in range(1, numNonZero + 1):
            colCoefList[ia[i] - 1] = da[i]
        return sparseList(colCoefList)
    
    def getRowCoef(self, rowIndex):
        """Returns a list of the row coefficients specified by rowIndex."""
        self._checkRowIndexValidity(rowIndex)
        ia = intArray(self.getNumCols() + 1)
        da = doubleArray(self.getNumCols() + 1)
        numNonZero = glp_get_mat_row(self.lp, rowIndex, ia, da)
        rowCoefList = [0. for i in range(self.getNumCols())]
        for i in range(1, numNonZero + 1):
            rowCoefList[ia[i] - 1] = da[i]
        return sparseList(rowCoefList)
    
    def changeColumn(self):
        """docstring for changeColumn"""
        # TODO not implemented
        pass
    
    def getObjective(self):
        """Returns a list of the objective coefficients."""
        num = self.getNumCols()
        return [glp_get_obj_coef(self.lp, i) for i in range(1, num + 1)]
    
    def setObjective(self, objDict):
        """objDict is a sparse dict of the objective function. The keys are the
        indices and the values the coefficients."""
        command = SetObjectiveCommand(self, objDict)
        command.execute()
    
    def getOptFlag(self):
        """Returns the optimization direction."""
        return glp_get_obj_dir(self.lp)
    
    def setOptFlag(self, optFlag='MAX'):
        """Sets optimization direction."""
        flag = ''.join([char.capitalize() for char in optFlag])
        command = SetOptFlagCommand(self, flag)
        command.execute()
    
    def primalValues(self):
        """Returns a list of all current primal values."""
        num = self.getNumCols()
        return [glp_get_col_prim(self.lp, i) for i in range(1, num + 1)]

    def dualValues(self):
        """Returns a list of all current primal values."""
        num = self.getNumCols()
        return [glp_get_col_dual(self.lp, i) for i in range(1, num + 1)]
    
    def getColumnsOfType(self, desiredKind):
        "Returns a list of column indices that have"
        if isinstance(desiredKind, str):
            mapping = {"integer":GLP_IV, "binary":GLP_BV, "continuous":GLP_CV}
            try:
                desiredKind = mapping[desiredKind]
            except KeyError, msg:
                raise KeyError, "Unknown column type " + desiredKind + '. Choose between "binary", "integer" and "continous".'
                sys.exit(-1)
        return [index for index, kind in self.getColumnKinds().items() if kind == desiredKind]
    
    def mipValues(self):
        mipIndices = self.getColumnsOfType(GLP_BV)
        mipIndices.extend(self.getColumnsOfType(GLP_IV))
        mipVals = dict()
        for index in mipIndices:
            mipVals[self.translateColumnIndices((index,))[0]] = glp_mip_col_val(self.lp, index)
        return mipVals
    
    def dualValues(self):
        """Returns a list of all current dual values."""
        num = self.getNumCols()
        return [glp_get_col_dual(self.lp, i) for i in range(1, num + 1)]
        
    def addPositiveColumnSwitches(self, columns=None):
        """
        If integer switch switch_c = 1 than lower_bound <= c <= upper_bound
        """
        if not columns:
            columns = self.getColumnIDs()
        cols = dict()
        for c in columns:
            cols["switch_"+c] = (0., 1., dict())
        # print 'adding columns'
        self.addColumns(cols)
        # print 'added columns'
        colKinds = dict()
        for c in columns:
            colKinds[self.translateColumnNames(("switch_"+c,))[0]] = GLP_BV
        self.setColumnKinds(colKinds)
        colBounds = self.getColumnBounds()
        newRows = dict()
        for c in columns:
            # print "processing ", r
            viPos = self.translateColumnNames((c,))[0]
            yPos = self.translateColumnNames(("switch_"+c,))[0]
            coefUb = {viPos:1., yPos:-colBounds[c][1]}
            coefLb = {viPos:1., yPos:-colBounds[c][0]}
            newRows["Ub_switch_"+c] = ('-inf', 0., coefUb)
            newRows["Lb_switch_"+c] = (0., 'inf', coefLb)
        # print 'adding rows'
        self.addRows(newRows)
        # print 'added rows'
    

    def addNegativeColumnSwitches(self, columns=None):
        """
        If integer switch switch_c = 1 than c = 0
        """
        if not columns:
            columns = self.getColumnIDs()
        cols = dict()
        for c in columns:
            cols["switch_"+c] = (0., 1., dict())
        cols["one_stub"] = (1., 1., dict())
        self.addColumns(cols)
        colKinds = dict()
        for c in columns:
            colKinds[self.translateColumnNames(("switch_"+c,))[0]] = GLP_BV
        self.setColumnKinds(colKinds)
        colBounds = self.getColumnBounds()
        newRows = dict()
        for c in columns:
            stubPos = self.translateColumnNames(("one_stub",))[0]
            viPos = self.translateColumnNames((c,))[0]
            yPos = self.translateColumnNames(("switch_"+c,))[0]
            coefUb = {viPos:1., yPos:colBounds[c][1], stubPos:-colBounds[c][1]}
            coefLb = {viPos:1., yPos:colBounds[c][0], stubPos:-colBounds[c][0]}
            newRows["Ub_switch_"+c] = ('-inf', 0., coefUb)
            newRows["Lb_switch_"+c] = (0., 'inf', coefLb)
        self.addRows(newRows)


class Command(object):
    """Command pattern stub. Defines only the interface."""
    def __init__(self, reciever):
        self.reciever = reciever
        self.memory = {}
    
    def execute(self):
        """docstring for execute"""
        pass
    
    def undo(self):
        """docstring for undo"""
        pass


class SetColumnBoundsCommand(Command):
    """Set bounds on reciever.
    
    Allows undoable modifications.
    """
    def __init__(self, reciever, boundDict):
        super(SetColumnBoundsCommand, self).__init__(reciever)
        self.boundDict = boundDict
        self.lp = self.reciever.lp
        for elem in self.boundDict:
            self.memory[elem] = (glp_get_col_lb(self.lp, elem),
                glp_get_col_ub(self.lp, elem))
    
    def execute(self):
        """docstring for execute"""
        self.reciever._modifyColumnBoundsNoMemory(self.boundDict)
        self.reciever.history.insert(0, self)
    
    def undo(self):
        self.reciever._modifyColumnBoundsNoMemory(self.memory)


class SetRowBoundsCommand(Command):
    """Set bounds on reciever.
    
    Allows undoable modifications.
    """
    def __init__(self, reciever, boundDict):
        super(SetRowBoundsCommand, self).__init__(reciever)
        self.boundDict = boundDict
        self.lp = self.reciever.lp
        for elem in self.boundDict:
            self.memory[elem] = (glp_get_row_lb(self.lp, elem),
                glp_get_row_ub(self.lp, elem))
    
    def execute(self):
        """docstring for execute"""
        self.reciever._modifyRowBoundsNoMemory(self.boundDict)
        self.reciever.history.insert(0, self)
    
    def undo(self):
        self.reciever._modifyRowBoundsNoMemory(self.memory)


class SetObjectiveCommand(Command):
    """Set new objective on reciever.
    
    objDict is a sparse dict of the objective function.
    The keys are the indices and the values the coefficients.
    Allows undoable modifications.
    """
    def __init__(self, reciever, objDict):
        super(SetObjectiveCommand, self).__init__(reciever)
        self.objDict = objDict
        self.memory = self.reciever.getObjective()
    
    def execute(self):
        """docstring for execute"""
        num = self.reciever.getNumCols()
        for i in range(1, num + 1):
            try:
                coeff = self.objDict[i]
                glp_set_obj_coef(self.reciever.lp, i, coeff)
            except KeyError:
                glp_set_obj_coef(self.reciever.lp, i, 0.)
        self.reciever.history.insert(0, self)
    
    def undo(self):
        for index, coef in enumerate(self.memory):
            glp_set_obj_coef(self.reciever.lp, index + 1, coef)


class SetOptFlagCommand(Command):
    """Set new optimization flag.
    
    Allows undoable modifications.
    """
    def __init__(self, reciever, flag):
        super(SetOptFlagCommand, self).__init__(reciever)
        self.flag = flag
        self.memory = self.reciever.getOptFlag()
    
    def execute(self):
        """docstring for execute"""
        glp_set_obj_dir(self.reciever.lp, eval("GLP_" + self.flag))
        self.reciever.history.insert(0, self)
    
    def undo(self):
        glp_set_obj_dir(self.reciever.lp, self.memory)


class DeleteColumnsCommand(Command):
    """Set new optimization flag.
    
    Allows undoable modifications. The order of the columns is not preserved
    """
    def __init__(self, reciever, colIndexes):
        super(DeleteColumnsCommand, self).__init__(reciever)
        self.colIndexes = colIndexes
        self.memory = self._memoryInit()
    
    def _memoryInit(self):
        """docstring for _memoryInit"""
        memory = dict()
        for elem in self.colIndexes:
            # get column id
            identifier = glp_get_col_name(self.reciever.lp, elem)
            # get lower and upper bounds
            lb = glp_get_col_lb(self.reciever.lp, elem)
            ub = glp_get_col_ub(self.reciever.lp, elem)
            colCoefList = self.reciever.getColumnCoef(elem)
            memory[identifier] = (lb, ub, colCoefList)
        return memory
    
    def execute(self):
        """docstring for execute"""
        l = len(self.colIndexes)
        num = intArray(l)
        for inc, i in enumerate(self.colIndexes):
            num[inc + 1] = i
        glp_del_cols(self.reciever.lp, l, num)
        self.reciever.history.insert(0, self)
    
    def undo(self):
        self.reciever.addColumns(self.memory)
        # The next line prevents that the undo step itself is recordes as an
        # undoable step
        self.reciever.history.pop(0)


class DeleteRowsCommand(Command):
    """Set new optimization flag.
    
    Allows undoable modifications. The order of the columns is not preserved
    """
    def __init__(self, reciever, rowsIndexes):
        super(DeleteRowsCommand, self).__init__(reciever)
        self.rowsIndexes = rowsIndexes
        self.memory = self._memoryInit()
    
    def _memoryInit(self):
        """docstring for _memoryInit"""
        memory = dict()
        for elem in self.rowsIndexes:
            # get row id
            identifier = glp_get_row_name(self.reciever.lp, elem)
            # get lower and upper bounds
            lb = glp_get_row_lb(self.reciever.lp, elem)
            ub = glp_get_row_ub(self.reciever.lp, elem)
            rowCoefList = self.reciever.getRowCoef(elem)
            memory[identifier] = (lb, ub, rowCoefList)
        return memory
    
    def execute(self):
        """docstring for execute"""
        l = len(self.rowsIndexes)
        num = intArray(l)
        for inc, i in enumerate(self.rowsIndexes):
           num[inc + 1] = i
        glp_del_rows(self.reciever.lp, l, num)
        self.reciever.history.insert(0, self)
    
    def undo(self):
        self.reciever.addRows(self.memory)
        # The next line prevents that the undo step itself is recordes as an
        # undoable step
        self.reciever.history.pop(0)


class AddColumnsCommand(Command):
    """Adds new columns to the lp.
    Allows undoable modifications.
    columns is dictionary of the form {columnID : (lb, ub, coeffList)}
    """
    def __init__(self, reciever, columns):
        super(AddColumnsCommand, self).__init__(reciever)
        self.columns = columns
        self.numCols = self.reciever.getNumCols()
        self.memory = range(self.numCols + 1, self.numCols+len(columns)+1)
    
    def execute(self):
        numCols = self.reciever.getNumCols()
        numRows = self.reciever.getNumRows()
        glp_add_cols(self.reciever.lp, len(self.columns))
        for i, item in enumerate(self.columns):
            index = numCols+1+i
            lb = self.columns[item][0]
            ub = self.columns[item][1]
            ia = intArray(numRows+1)
            da = doubleArray(numRows+1)
            sparseCoef = self.columns[item][2]
            for j in range(1, numRows + 1):
                if j in sparseCoef:
                    ia[j] = j
                    da[j] = sparseCoef[j]
                else:
                    ia[j] = j
                    da[j] = 0.
            glp_set_col_name(self.reciever.lp, index, item)
            glp_set_obj_coef(self.reciever.lp, index, 0.)
            self.reciever._setColumnBound(index, lb, ub)
            glp_set_mat_col(self.reciever.lp, index, numRows, ia, da)
        self.reciever.history.insert(0, self)
    
    def undo(self):
        self.reciever.deleteColumns(self.memory)
        # The next line prevents that the undo step itself is recordes as an
        # undoable step
        self.reciever.history.pop(0)


class SetColumnKindsCommand(Command):
    """Sets the kinds of columns
    Allows undoable modifications.
    kindsDict is dictionary of the form {columnIndex : columnKind, ...}
    """
    def __init__(self, reciever, kindsDict):
        super(SetColumnKindsCommand, self).__init__(reciever)
        self.kindsDict = kindsDict
        self.memory = (self.reciever.getColumnKinds(kindsDict.keys()), self.reciever.getColumnBounds())
    
    def execute(self):
        for i in self.kindsDict:
            glp_set_col_kind(self.reciever.lp, i, self.kindsDict[i])
        self.reciever.history.insert(0, self)
    
    def undo(self):
        self.reciever.setColumnKinds(self.memory[0])
        # The next line prevents that the undo step itself is recorded as an
        # undoable step
        self.reciever.history.pop(0)
        self.reciever.modifyColumnBounds(self.memory[1])
        self.reciever.history.pop(0)


class AddRowsCommand(Command):
    """Adds new rows to the lp.
    Allows undoable modifications.
    rows is dictionary of the form {rowID : (lb, ub, coeffList)}
    """
    def __init__(self, reciever, rows):
        super(AddRowsCommand, self).__init__(reciever)
        self.rows = rows
        self.numRows = self.reciever.getNumRows()
        self.memory = range(self.numRows + 1, self.numRows+len(rows)+1)
    
    def execute(self):
        numCols = self.reciever.getNumCols()
        numRows = self.reciever.getNumRows()
        glp_add_rows(self.reciever.lp, len(self.rows))
        for i, item in enumerate(self.rows):
            index = numRows+1+i
            lb = self.rows[item][0]
            ub = self.rows[item][1]
            ia = intArray(numCols+1)
            da = doubleArray(numCols+1)
            sparseCoef = self.rows[item][2]
            for j in range(1, numCols + 1):
                if j in sparseCoef:
                    ia[j] = j
                    da[j] = sparseCoef[j]
                else:
                    ia[j] = j
                    da[j] = 0.
            glp_set_row_name(self.reciever.lp, index, item)
            # glp_set_obj_coef(self.reciever.lp, index, 0.)
            self.reciever._setRowBound(index, lb, ub)
            glp_set_mat_row(self.reciever.lp, index, numCols, ia, da)
        self.reciever.history.insert(0, self)
    
    def undo(self):
        self.reciever.deleteRows(self.memory)
        # The next line prevents that the undo step itself is recordes as an
        # undoable step
        self.reciever.history.pop(0)


class sparseList(dict):
    
    def __init__(self, l):
        """docstring for __init__"""
        super(sparseList, self).__init__(self._initList(l))
    
    def _initList(self, l):
        tmp = list()
        for index, val in enumerate(l):
            if val != 0:
                tmp.append((index + 1, val))
        return tmp
    


if __name__ == '__main__':
    import util
    import pprint
    lp = glpk(util.ImportCplex('test_data/model.lp'))
    # lp.toggleVerbosity()
    # print lp.getColumnKinds((1,))
    #     lp.setColumnKinds({1:GLP_BV})
    #     print lp.getColumnKinds((1,))
    #     lp.initialize()
    #     print lp.getColumnKinds((1,))
    #     print lp.getColumnBounds()
    #     lp.setColumnKinds({1:1})
    #     print lp.getColumnKinds((1,))
    #     lp.initialize()
    #     print lp.getColumnKinds((1,))
    lp.addRows({'R("R_HansWurs")': (0., 99999., {})})
    lp.simplex()
    print lp.getObjVal()