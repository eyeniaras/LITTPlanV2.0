/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was created and modified by Erol Yeniaras using the template originally 
  developed by Jean-Christophe Fillion-Robin.

==============================================================================*/
#include "ctkUtils.h"

// Qt includes
#include <QFileDialog>

// SlicerQt includes
#include "qSlicerLITTPlanV2ModuleWidget.h"
#include "ui_qSlicerLITTPlanV2Module.h"
//#include "qSlicerApplication.h"
//#include "qSlicerIOManager.h"

// vtkSlicerLogic includes
#include "vtkSlicerTransformLogic.h"

// MRMLWidgets includes
#include <qMRMLUtils.h>

// MRML includes
#include "vtkMRMLLinearTransformNode.h"
#include "vtkMRMLAnnotationHierarchyNode.h"
#include "vtkMRMLAnnotationFiducialNode.h"
#include "vtkCollection.h"
#include "vtkMRMLFiducial.h"
#include "vtkMRMLFiducialListNode.h"
#include "vtkMRMLFiducialListStorageNode.h"
#include "vtkMRMLLinearTransformNode.h"
#include "vtkMRMLScene.h"
#include "vtkMRMLDisplayNode.h"
#include "vtkMRMLModelDisplayNode.h"


// VTK includes
#include "vtkPolyData.h"
#include <vtkLineSource.h>
#include <vtkSphereSource.h>
#include <vtkCylinderSource.h>
#include <vtkCellArray.h>

#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>

#include "vtkSTLReader.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataNormals.h"
#include <vtkExodusIIReader.h>

//-----------------------------------------------------------------------------
class qSlicerLITTPlanV2ModuleWidgetPrivate: public Ui_qSlicerLITTPlanV2Module
{
  Q_DECLARE_PUBLIC(qSlicerLITTPlanV2ModuleWidget);
protected:
  qSlicerLITTPlanV2ModuleWidget* const q_ptr;
public:
  qSlicerLITTPlanV2ModuleWidgetPrivate(qSlicerLITTPlanV2ModuleWidget& object);
  vtkSlicerTransformLogic*      logic()const;
  QButtonGroup*                 CoordinateReferenceButtonGroup;
  vtkMRMLLinearTransformNode*   MRMLTransformNode;
};

static double* coors0;
static double* coors1;

//-----------------------------------------------------------------------------
qSlicerLITTPlanV2ModuleWidgetPrivate::qSlicerLITTPlanV2ModuleWidgetPrivate(qSlicerLITTPlanV2ModuleWidget& object)
  : q_ptr(&object)
{
  this->CoordinateReferenceButtonGroup = 0;
  this->MRMLTransformNode = 0;
}
//-----------------------------------------------------------------------------
vtkSlicerTransformLogic* qSlicerLITTPlanV2ModuleWidgetPrivate::logic()const
{
  Q_Q(const qSlicerLITTPlanV2ModuleWidget);
  return vtkSlicerTransformLogic::SafeDownCast(q->logic());
}

//-----------------------------------------------------------------------------
qSlicerLITTPlanV2ModuleWidget::qSlicerLITTPlanV2ModuleWidget(QWidget* _parentWidget)
  : Superclass(_parentWidget)
  , d_ptr(new qSlicerLITTPlanV2ModuleWidgetPrivate(*this))
{
}

//-----------------------------------------------------------------------------
qSlicerLITTPlanV2ModuleWidget::~qSlicerLITTPlanV2ModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::setup()
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  d->setupUi(this);

  // Add coordinate reference button to a button group
  d->CoordinateReferenceButtonGroup =
    new QButtonGroup(d->CoordinateReferenceGroupBox);
  d->CoordinateReferenceButtonGroup->addButton(
    d->GlobalRadioButton, qMRMLTransformSliders::GLOBAL);
  d->CoordinateReferenceButtonGroup->addButton(
    d->LocalRadioButton, qMRMLTransformSliders::LOCAL);

  connect(this, SIGNAL(mrmlSceneChanged(vtkMRMLScene*)), d->inputFiducialsNodeSelector, SLOT(setMRMLScene(vtkMRMLScene*)));
  connect(d->btnPath, SIGNAL(clicked()), this, SLOT(CreatePath()));
  connect(d->btnPoints, SIGNAL(clicked()), this, SLOT(ExportPoints()));
  connect(d->btnApplicator, SIGNAL(clicked()), this, SLOT(LoadApplicator()));
  connect(d->btnFEM, SIGNAL(clicked()), this, SLOT(GetFEM()));
  connect(d->btnBurn, SIGNAL(clicked()), this, SLOT(onBtnBurnClicked()));
  
  // Connect button group
  this->connect(d->CoordinateReferenceButtonGroup,
                SIGNAL(buttonPressed(int)),
                SLOT(onCoordinateReferenceButtonPressed(int)));

  // Connect identity button
  this->connect(d->IdentityPushButton,
                SIGNAL(clicked()),
                SLOT(identity()));

  // Connect revert button
  this->connect(d->InvertPushButton,
                SIGNAL(clicked()),
                SLOT(invert()));

  // Connect node selector with module itself
  this->connect(d->TransformNodeSelector,
                SIGNAL(currentNodeChanged(vtkMRMLNode*)),
                SLOT(onNodeSelected(vtkMRMLNode*)));

  // Connect minimum and maximum from the translation sliders to the matrix
  this->connect(d->TranslationSliders,
               SIGNAL(rangeChanged(double,double)),
               SLOT(onTranslationRangeChanged(double,double)));

  // Notify the matrix of the current translation min/max values
  this->onTranslationRangeChanged(d->TranslationSliders->minimum(),
                                  d->TranslationSliders->maximum());

  // Transform nodes connection
  this->connect(d->TransformToolButton, SIGNAL(clicked()),
                SLOT(transformSelectedNodes()));
  this->connect(d->UntransformToolButton, SIGNAL(clicked()),
                SLOT(untransformSelectedNodes()));

  // Icons
  QIcon rightIcon =
    QApplication::style()->standardIcon(QStyle::SP_ArrowRight);
  d->TransformToolButton->setIcon(rightIcon);

  QIcon leftIcon =
    QApplication::style()->standardIcon(QStyle::SP_ArrowLeft);
  d->UntransformToolButton->setIcon(leftIcon);

  this->onNodeSelected(0);
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::LoadApplicator()
{	
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("C:\\Users\\eyeniaras\\Downloads\\daviddata\\laserApplicator.stl");
	
	//vtkSmartPointer<vtkExodusIIReader> reader = vtkSmartPointer<vtkExodusIIReader>::New();
	//   reader->SetFileName("C:\\Users\\eyeniaras\\Downloads\\daviddata\\fem.e");
	//reader->ExodusModelMetadataOn();
	//reader->UpdateInformation();
	////reader->SetPointResultArrayStatus("u", 0);
	//vtkSmartPointer<vtkContourFilter> filter = vtkSmartPointer<vtkContourFilter>::New();
	//filter->SetInputConnection(reader->GetOutputPort());
	////filter->SetValue(0, 50);
	//vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
	//normals->SetInputConnection(filter->GetOutputPort());
	////normals->SetFeatureAngle(60.0);

	vtkMRMLLinearTransformNode* transformNode=vtkMRMLLinearTransformNode::New();
	this->mrmlScene()->AddNode(transformNode);

	vtkSmartPointer<vtkMatrix4x4> mat = vtkSmartPointer<vtkMatrix4x4>::New();
    mat->Identity();
	
	//mat->SetElement(0,3,coors0[0]);
	//mat->SetElement(1,3,coors0[1]);
	//mat->SetElement(2,3,coors0[2]);
    /*for( size_t p = 0; p < 4; p++ )
      {
      vtkFloatingPointType point = normal[p];
      mat->SetElement(static_cast<int>(p), 0, (sign * point) );
      }*/
    /*vtkFloatingPointType oneAndAlpha = 1.0 + mat->GetElement(0, 0);
    mat->SetElement(0, 1, -1.0 * mat->GetElement(1, 0) );
    mat->SetElement(0, 2, (-1.0 * (mat->GetElement(2, 0) ) ) );
    mat->SetElement(2, 1, (-1.0 * (mat->GetElement(1, 0) * (mat->GetElement(2, 0) / oneAndAlpha) ) ) );
    mat->SetElement(1, 2, (-1.0 * (mat->GetElement(1, 0) * (mat->GetElement(2, 0) / oneAndAlpha) ) ) );
    mat->SetElement(1, 1, (1.0  - (mat->GetElement(1, 0) * (mat->GetElement(1, 0) / oneAndAlpha) ) ) );
    mat->SetElement(2, 2, (1.0  - (mat->GetElement(2, 0) * (mat->GetElement(2, 0) / oneAndAlpha) ) ) );*/


	transformNode->SetAndObserveMatrixTransformToParent( mat.GetPointer() );

	vtkMRMLModelNode* modelNode=vtkMRMLModelNode::New();
	modelNode->SetAndObserveTransformNodeID(transformNode->GetID());
	modelNode->SetScene(this->mrmlScene());
	modelNode->SetName("Applicator");
	modelNode->SetAndObservePolyData(reader->GetOutput());	
    reader->Update();

	vtkMRMLModelDisplayNode* modelDisplayNode=vtkMRMLModelDisplayNode::New();
	modelDisplayNode->SetColor(1,0,0); // green
	modelDisplayNode->SetScene(this->mrmlScene());
	this->mrmlScene()->AddNode(modelDisplayNode);
	modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

	modelDisplayNode->SetPolyData(modelNode->GetPolyData());
	this->mrmlScene()->AddNode(modelNode);

	transformNode->Delete();
	modelNode->Delete();
	modelDisplayNode->Delete();	    
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::GetFEM()
{
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::ExportPoints()
{
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::CreatePath()
{	
	// Alternate Approach
	//vtkMRMLNode* fiducialsNode = d_ptr->inputFiducialsNodeSelector->currentNode();

	
	//vtkMRMLFiducialListNode *fnode=vtkMRMLFiducialListNode::SafeDownCast(fiducialsNode); 
	
	//vtkMRMLAnnotationHierarchyNode *anode=vtkMRMLAnnotationHierarchyNode::SafeDownCast(fiducialsNode);

	//vtkCollection *cnodes = vtkCollection::New();	

	////if(anode->GetClassName()== "vtkMRMLAnnotationHierarchyNode")
	////{
	////	//d_ptr->label->clear();
	//	anode->GetChildrenDisplayableNodes(cnodes);

	////	
	//	int n=cnodes->GetNumberOfItems();
	//	if (n >0)
	//	{

	//		for(int i=0; i<n; i++)
	//		{
	//			d_ptr->pushButton->setText(QString::number(n));
	//			

	//			vtkMRMLFiducial *fid = vtkMRMLFiducial::SafeDownCast(cnodes->GetItemAsObject(i));
	//			if (fid != NULL)
	//			{
	//				//d_ptr->label->clear();
	//				float* xyz = fid->GetXYZ();
	//				//d_ptr->pushButton->setText(QString::number(xyz[0]));
	//			}
	//			//float coords[];
	//			//fid->GetXYZ(coords);
	//			
	//			/*vtkMRMLFiducial *fid2 = vtkMRMLFiducial::SafeDownCast(cnodes->GetItemAsObject(1));
	//			float coords2[3];
	//			fid1->GetXYZ(coords2);

	//			
	//			QString qnum= QString::number(n);
	//			d_ptr->pushButton->setText(QString::number(coords1[0]));*/
	//		}
	//	}
	////}
	//cnodes->RemoveAllItems();
	//cnodes->Delete();
	////
	// // OR Get nodes from the mrml scene...
	
	std::vector<vtkMRMLNode*> nodes;
	this->mrmlScene()->GetNodesByClass("vtkMRMLAnnotationHierarchyNode", nodes);
	
	if (nodes.size() > 1 && (d_ptr->inputFiducialsNodeSelector->currentNode() != 0))
	{		
		for (unsigned int i=0; i<nodes.size(); i++)
		{
			vtkMRMLAnnotationHierarchyNode *hnode = vtkMRMLAnnotationHierarchyNode::SafeDownCast(nodes[i]);
			vtkCollection *cnodes = vtkCollection::New();
			hnode->GetDirectChildren(cnodes);
			if (cnodes->GetNumberOfItems() > 0 && vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0)) != NULL)
			{
				vtkMRMLAnnotationFiducialNode *fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0));
				coors0=fnode0->GetFiducialCoordinates();
				vtkMRMLAnnotationFiducialNode *fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1));
				coors1=fnode1->GetFiducialCoordinates();
								
				vtkLineSource* line=vtkLineSource::New();
				line->SetPoint1(coors0[0], coors0[1], coors0[2]);
				line->SetPoint2(coors1[0], coors1[1], coors1[2]);
							
				vtkPolyData* polydata = vtkPolyData::New();				
				polydata=line->GetOutput();
				
				vtkMRMLModelNode* modelNode=vtkMRMLModelNode::New();
				modelNode->SetScene(this->mrmlScene());
				modelNode->SetName("ApplicatorPath");
				modelNode->SetAndObservePolyData(polydata);

				vtkMRMLModelDisplayNode* modelDisplayNode=vtkMRMLModelDisplayNode::New();
				modelDisplayNode->SetColor(0,1,0); // green
				modelDisplayNode->SetScene(this->mrmlScene());
				this->mrmlScene()->AddNode(modelDisplayNode);
				modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

				modelDisplayNode->SetPolyData(modelNode->GetPolyData());
				this->mrmlScene()->AddNode(modelNode);
												
				line->Delete();
				modelNode->Delete();
				modelDisplayNode->Delete();
				polydata->Delete();
				cnodes->RemoveAllItems();
				cnodes->Delete();
				break;
			}
			cnodes->RemoveAllItems();
			cnodes->Delete();
		}
	}
	nodes.clear();	
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onBtnApplicatorClicked()
{
	d_ptr->btnApplicator->setText("alo");
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("C:\\Users\\eyeniaras\\Downloads\\daviddata\\laserApplicator.stl");
	
	//vtkSmartPointer<vtkExodusIIReader> reader = vtkSmartPointer<vtkExodusIIReader>::New();
 //   reader->SetFileName("C:\\Users\\eyeniaras\\Downloads\\daviddata\\fem.e");
	//reader->ExodusModelMetadataOn();
	//reader->UpdateInformation();
	////reader->SetPointResultArrayStatus("u", 0);
	//vtkSmartPointer<vtkContourFilter> filter = vtkSmartPointer<vtkContourFilter>::New();
	//filter->SetInputConnection(reader->GetOutputPort());
	////filter->SetValue(0, 50);
	//vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
	//normals->SetInputConnection(filter->GetOutputPort());
	////normals->SetFeatureAngle(60.0);

	vtkMRMLLinearTransformNode* transformNode=vtkMRMLLinearTransformNode::New();
	this->mrmlScene()->AddNode(transformNode);

	vtkSmartPointer<vtkMatrix4x4> mat = vtkSmartPointer<vtkMatrix4x4>::New();
    mat->Identity();
	//mat->SetElement(2,3,86.0);
    /*for( size_t p = 0; p < 4; p++ )
      {
      vtkFloatingPointType point = normal[p];
      mat->SetElement(static_cast<int>(p), 0, (sign * point) );
      }*/
    /*vtkFloatingPointType oneAndAlpha = 1.0 + mat->GetElement(0, 0);
    mat->SetElement(0, 1, -1.0 * mat->GetElement(1, 0) );
    mat->SetElement(0, 2, (-1.0 * (mat->GetElement(2, 0) ) ) );
    mat->SetElement(2, 1, (-1.0 * (mat->GetElement(1, 0) * (mat->GetElement(2, 0) / oneAndAlpha) ) ) );
    mat->SetElement(1, 2, (-1.0 * (mat->GetElement(1, 0) * (mat->GetElement(2, 0) / oneAndAlpha) ) ) );
    mat->SetElement(1, 1, (1.0  - (mat->GetElement(1, 0) * (mat->GetElement(1, 0) / oneAndAlpha) ) ) );
    mat->SetElement(2, 2, (1.0  - (mat->GetElement(2, 0) * (mat->GetElement(2, 0) / oneAndAlpha) ) ) );*/


	transformNode->SetAndObserveMatrixTransformToParent( mat.GetPointer() );

	vtkMRMLModelNode* modelNode=vtkMRMLModelNode::New();
	modelNode->SetAndObserveTransformNodeID(transformNode->GetID());
	modelNode->SetScene(this->mrmlScene());
	modelNode->SetName("Applicator");
	modelNode->SetAndObservePolyData(reader->GetOutput());	
    reader->Update();

	vtkMRMLModelDisplayNode* modelDisplayNode=vtkMRMLModelDisplayNode::New();
	modelDisplayNode->SetColor(1,0,0); // green
	modelDisplayNode->SetScene(this->mrmlScene());
	this->mrmlScene()->AddNode(modelDisplayNode);
	modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

	modelDisplayNode->SetPolyData(modelNode->GetPolyData());
	this->mrmlScene()->AddNode(modelNode);

	transformNode->Delete();
	modelNode->Delete();
	modelDisplayNode->Delete();	    
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onBtnBurnClicked()
{
	vtkSphereSource* sphere=vtkSphereSource::New();
	sphere->SetRadius(20.0);
	sphere->SetCenter(.0, .0,.0);
							
	vtkPolyData* polydata = vtkPolyData::New();				
	polydata=sphere->GetOutput();
				
	vtkMRMLModelNode* modelNode=vtkMRMLModelNode::New();
	modelNode->SetScene(this->mrmlScene());
	modelNode->SetName("ApplicatorPath");
	modelNode->SetAndObservePolyData(polydata);

	vtkMRMLModelDisplayNode* modelDisplayNode=vtkMRMLModelDisplayNode::New();
	modelDisplayNode->SetColor(0,1,0); // green
	modelDisplayNode->SetScene(this->mrmlScene());
	this->mrmlScene()->AddNode(modelDisplayNode);
	modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

	modelDisplayNode->SetPolyData(modelNode->GetPolyData());
	this->mrmlScene()->AddNode(modelNode);
											
	sphere->Delete();
	modelNode->Delete();
	modelDisplayNode->Delete();
	polydata->Delete();
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onCoordinateReferenceButtonPressed(int id)
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  
  qMRMLTransformSliders::CoordinateReferenceType ref =
    (id == qMRMLTransformSliders::GLOBAL) ? qMRMLTransformSliders::GLOBAL : qMRMLTransformSliders::LOCAL;
  d->TranslationSliders->setCoordinateReference(ref);
  d->RotationSliders->setCoordinateReference(ref);
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onNodeSelected(vtkMRMLNode* node)
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  
  vtkMRMLLinearTransformNode* transformNode = vtkMRMLLinearTransformNode::SafeDownCast(node);

  // Enable/Disable CoordinateReference, identity buttons, MatrixViewGroupBox,
  // Min/Max translation inputs
  d->CoordinateReferenceGroupBox->setEnabled(transformNode != 0);
  d->IdentityPushButton->setEnabled(transformNode != 0);
  d->InvertPushButton->setEnabled(transformNode != 0);
  d->MatrixViewGroupBox->setEnabled(transformNode != 0);

  // Listen for Transform node changes
  this->qvtkReconnect(d->MRMLTransformNode, transformNode,
    vtkMRMLTransformableNode::TransformModifiedEvent,
    this, SLOT(onMRMLTransformNodeModified(vtkObject*)));

  QStringList nodeTypes;
  // If no transform node, it would show the entire scene, lets shown none
  // instead.
  if (transformNode == 0)
    {
    nodeTypes << QString("vtkMRMLNotANode");
    }
  d->TransformedTreeView->setNodeTypes(nodeTypes);

  // Filter the current node in the transformed tree view
  d->TransformedTreeView->setRootNode(transformNode);

  // Hide the current node in the transformable tree view
  QStringList hiddenNodeIDs;
  if (transformNode)
    {
    hiddenNodeIDs << QString(transformNode->GetID());
    }
  d->TransformableTreeView->sortFilterProxyModel()
    ->setHiddenNodeIDs(hiddenNodeIDs);
  d->MRMLTransformNode = transformNode;
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::identity()
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);

  if (!d->MRMLTransformNode)
    {
    return;
    }

  d->RotationSliders->resetUnactiveSliders();
  d->MRMLTransformNode->GetMatrixTransformToParent()->Identity();
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::invert()
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  
  if (!d->MRMLTransformNode) { return; }

  d->RotationSliders->resetUnactiveSliders();
  d->MRMLTransformNode->GetMatrixTransformToParent()->Invert();
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onMRMLTransformNodeModified(vtkObject* caller)
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  
  vtkMRMLLinearTransformNode* transformNode = vtkMRMLLinearTransformNode::SafeDownCast(caller);
  if (!transformNode) { return; }

  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  qMRMLUtils::getTransformInCoordinateSystem(d->MRMLTransformNode,
    this->coordinateReference() == qMRMLTransformSliders::GLOBAL, transform);

  // The matrix can be changed externally. The min/max values shall be updated 
  //accordingly to the new matrix if needed.
  vtkMatrix4x4 * mat = transform->GetMatrix();
  double min = 0.;
  double max = 0.;
  this->extractMinMaxTranslationValue(mat, min, max);
  if (min < d->TranslationSliders->minimum())
    {
    min = min - 0.3 * fabs(min);
    d->TranslationSliders->setMinimum(min);
    }
  if (max > d->TranslationSliders->maximum())
    {
    max = max + 0.3 * fabs(max);
    d->TranslationSliders->setMaximum(max);
    }
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::extractMinMaxTranslationValue(
  vtkMatrix4x4 * mat, double& min, double& max)
{
  if (!mat)
    {
    Q_ASSERT(mat);
    return;
    }
  for (int i=0; i <3; i++)
    {
    min = qMin(min, mat->GetElement(i,3));
    max = qMax(max, mat->GetElement(i,3));
    }
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onTranslationRangeChanged(double newMin,
                                                              double newMax)
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  d->MatrixWidget->setRange(newMin, newMax);
}

//-----------------------------------------------------------------------------
int qSlicerLITTPlanV2ModuleWidget::coordinateReference()const
{
  Q_D(const qSlicerLITTPlanV2ModuleWidget);
  return d->CoordinateReferenceButtonGroup->checkedId();
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::setMRMLScene(vtkMRMLScene* scene)
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  this->Superclass::setMRMLScene(scene);
  // If the root index is set before the scene, it will show the scene as
  // top-level item. Setting the root index to be the scene makes the nodes
  // top-level, and this can only be done once the scene is set.
  d->TransformableTreeView->setRootIndex(
    d->TransformableTreeView->sortFilterProxyModel()->mrmlSceneIndex());
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::transformSelectedNodes()
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  QModelIndexList selectedIndexes =
    d->TransformableTreeView->selectionModel()->selectedRows();
  selectedIndexes = qMRMLTreeView::removeChildren(selectedIndexes);
  foreach(QModelIndex selectedIndex, selectedIndexes)
    {
    vtkMRMLTransformableNode* node = vtkMRMLTransformableNode::SafeDownCast(
    d->TransformableTreeView->sortFilterProxyModel()->
      mrmlNodeFromIndex( selectedIndex ));
    Q_ASSERT(node);
    node->SetAndObserveTransformNodeID(d->MRMLTransformNode->GetID());
    }
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::untransformSelectedNodes()
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  QModelIndexList selectedIndexes =
    d->TransformedTreeView->selectionModel()->selectedRows();
  selectedIndexes = qMRMLTreeView::removeChildren(selectedIndexes);
  foreach(QModelIndex selectedIndex, selectedIndexes)
    {
    vtkMRMLTransformableNode* node = vtkMRMLTransformableNode::SafeDownCast(
    d->TransformedTreeView->sortFilterProxyModel()->
      mrmlNodeFromIndex( selectedIndex ));
    Q_ASSERT(node);
    node->SetAndObserveTransformNodeID(0);
    }
}
