#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PV4FoamReader'
viscoelasticOpenFOAM = PV4FoamReader(FileName='/home/mike/foam/mike-4.0/run/v2f/Channel/Viscoelastic/Viscoelastic.OpenFOAM')
viscoelasticOpenFOAM.MeshParts = ['internalMesh']
viscoelasticOpenFOAM.VolumeFields = ['p', 'U', 'C', 'RStress', 'yPlus', 'NLTPlus']

animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
viscoelasticOpenFOAMDisplay = Show(viscoelasticOpenFOAM, renderView1)
viscoelasticOpenFOAMDisplay.ColorArrayName = [None, '']
viscoelasticOpenFOAMDisplay.ScalarOpacityUnitDistance = 0.19007996047274922

# reset view to fit data
renderView1.ResetCamera()

# set scalar coloring
ColorBy(viscoelasticOpenFOAMDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
viscoelasticOpenFOAMDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# create a new 'Calculator'
calculator1 = Calculator(Input=viscoelasticOpenFOAM)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'sqrtRxx'
calculator1.Function = 'sqrt(volPointInterpolate(RStress)_XX)'
sqrtRxxLUT = GetColorTransferFunction('sqrtRxx')

# show data in view
calculator1Display = Show(calculator1, renderView1)
# trace defaults for the display properties.
calculator1Display.ColorArrayName = ['POINTS', 'sqrtRxx']
calculator1Display.LookupTable = sqrtRxxLUT
calculator1Display.ScalarOpacityUnitDistance = 0.19007996047274922

# hide data in view
Hide(viscoelasticOpenFOAM, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'sqrtRxx'
sqrtRxxPWF = GetOpacityTransferFunction('sqrtRxx')

# create a new 'Calculator'
calculator2 = Calculator(Input=calculator1)
calculator2.Function = ''

# Properties modified on calculator2
calculator2.ResultArrayName = 'sqrtRyy'
calculator2.Function = 'sqrt(volPointInterpolate(RStress)_YY)'

# get color transfer function/color map for 'sqrtRyy'
sqrtRyyLUT = GetColorTransferFunction('sqrtRyy')

# show data in view
calculator2Display = Show(calculator2, renderView1)
# trace defaults for the display properties.
calculator2Display.ColorArrayName = ['POINTS', 'sqrtRyy']
calculator2Display.LookupTable = sqrtRyyLUT
calculator2Display.ScalarOpacityUnitDistance = 0.19007996047274922

# hide data in view
Hide(calculator1, renderView1)

# show color bar/color legend
calculator2Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'sqrtRyy'
sqrtRyyPWF = GetOpacityTransferFunction('sqrtRyy')

# create a new 'Calculator'
calculator3 = Calculator(Input=calculator2)
calculator3.Function = ''

# Properties modified on calculator3
calculator3.ResultArrayName = 'sqrtRzz'
calculator3.Function = 'sqrt(volPointInterpolate(RStress)_ZZ)'

# get color transfer function/color map for 'sqrtRzz'
sqrtRzzLUT = GetColorTransferFunction('sqrtRzz')

# show data in view
calculator3Display = Show(calculator3, renderView1)
# trace defaults for the display properties.
calculator3Display.ColorArrayName = ['POINTS', 'sqrtRzz']
calculator3Display.LookupTable = sqrtRzzLUT
calculator3Display.ScalarOpacityUnitDistance = 0.19007996047274922

# hide data in view
Hide(calculator2, renderView1)

# show color bar/color legend
calculator3Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'sqrtRzz'
sqrtRzzPWF = GetOpacityTransferFunction('sqrtRzz')

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(Input=calculator3,
    Source='High Resolution Line Source')

# init the 'High Resolution Line Source' selected for 'Source'
plotOverLine1.Source.Point2 = [0.10000000149011612, 1.0, 0.10000000149011612]

# Properties modified on plotOverLine1
plotOverLine1.Tolerance = 2.22044604925031e-16

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Point1 = [0.05000000074505806, 0.0, 0.05000000074505806]
plotOverLine1.Source.Point2 = [0.05000000074505806, 1.0, 0.05000000074505806]

# get layout
viewLayout1 = GetLayout()

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')
lineChartView1.ViewSize = [658, 860]

# place view in the layout
viewLayout1.AssignView(2, lineChartView1)

# show data in view
plotOverLine1Display = Show(plotOverLine1, lineChartView1)
# trace defaults for the display properties.
plotOverLine1Display.CompositeDataSetIndex = [0]
plotOverLine1Display.UseIndexForXAxis = 0
plotOverLine1Display.XArrayName = 'arc_length'

# save data
SaveData('/home/mike/foam/mike-4.0/run/v2f/Channel/Viscoelastic/DATA.csv', proxy=plotOverLine1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.05000000074505806, 0.5, 2.001074531854705]
renderView1.CameraFocalPoint = [0.05000000074505806, 0.5, 0.05000000074505806]
renderView1.CameraParallelScale = 0.5049752470656473

import sys
sys.exit()

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
