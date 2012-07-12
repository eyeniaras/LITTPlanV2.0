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
#include <fstream>
#include <time.h>
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
#include <vtkCallbackCommand.h>
#include "vtkPolyData.h"
#include <vtkLineSource.h>
#include <vtkSphereSource.h>
#include <vtkCylinderSource.h>
#include <vtkCellArray.h>

#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

#include "vtkSTLReader.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataNormals.h"
#include <vtkExodusIIReader.h>
//-----------------------------------------------------------------------------
  static int count=0;

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
	SetParametersToDefault();	
	tissueIndex=0;

  this->CallBack1 = vtkSmartPointer<vtkCallbackCommand>::New();
  this->CallBack2 = vtkSmartPointer<vtkCallbackCommand>::New();

  this->fnode0 = vtkSmartPointer<vtkMRMLAnnotationFiducialNode>::New();
  this->fnode1 = vtkSmartPointer<vtkMRMLAnnotationFiducialNode>::New();

  observersAreActive=false;
  fiducialsAreValid=false;
  applicatorVisible=true;
  pathVisible=true;
  startPoint[0]=0.0; startPoint[1]=0.0; startPoint[2]=0.0; 
  targetPoint[0]=0.0; targetPoint[1]=0.0; targetPoint[2]=0.0;
  prevStartPoint[0]=0.0; prevStartPoint[1]=0.0; prevStartPoint[2]=0.0; 
  prevTargetPoint[0]=0.0; prevTargetPoint[1]=0.0; prevTargetPoint[2]=0.0;
}

//-----------------------------------------------------------------------------
qSlicerLITTPlanV2ModuleWidget::~qSlicerLITTPlanV2ModuleWidget()
{	
	if (this->mrmlScene())
    {
		this->mrmlScene()->RemoveAllObservers();
		//this->mrmlScene()->RemoveObserver(this->CallBack1);
		//this->mrmlScene()->RemoveObserver(this->CallBack2);
    }
	this->observersAreActive=false;	
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::setup()
{
  Q_D(qSlicerLITTPlanV2ModuleWidget);
  d->setupUi(this);

  CallBack1->SetClientData(this);
  CallBack2->SetClientData(this);
  
  // Initialize the parameters textboxes
  d->ThermalConductivity->setText("0.00");
  d->TissuePerfusion->setText("0.00");
  d->OpticalAbsorption->setText("0.00");
  d->OpticalScattering->setText("0.00");
  d->OpticalAnisotrophy->setText("0.00");
  d_ptr->ThermalDose->setText(QString::number(thermalDose));
  d_ptr->Power->setText(QString::number(power));
 
  // When MRML scene changed load its contens to inputFiducials node selector (the list of the fiducials)
  connect(this, SIGNAL(mrmlSceneChanged(vtkMRMLScene*)), d->inputFiducialsNodeSelector, SLOT(setMRMLScene(vtkMRMLScene*)));

  // When the fiducials are selected load them into the class
  connect(d->inputFiducialsNodeSelector, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(onFiducalNodeChanged(vtkMRMLNode*)));
 
  // Hide Path button is clicked
  connect(d->btnHide, SIGNAL(clicked()), this, SLOT(HidePath()));

  // Show Path button is clicked
  connect(d->btnPath, SIGNAL(clicked()), this, SLOT(CreatePath()));
  
  // Load Applicator button is clicked
  connect(d->btnApplicator, SIGNAL(clicked()), this, SLOT(LoadApplicator()));
  
  //connect(d->btnFEM, SIGNAL(clicked()), this, SLOT(GetFEM()));
  connect(d->btnBurn, SIGNAL(clicked()), this, SLOT(onBtnBurnClicked()));
  connect(d->TissueSelector, SIGNAL(currentIndexChanged(int)), this, SLOT(onTissueTypeChanged()));
  connect(d->btnDefault, SIGNAL(clicked()), this, SLOT(onBtnDefaultClicked()));
  
  // Add coordinate reference button to a button group
  d->CoordinateReferenceButtonGroup = new QButtonGroup(d->CoordinateReferenceGroupBox);
  d->CoordinateReferenceButtonGroup->addButton(d->GlobalRadioButton, qMRMLTransformSliders::GLOBAL);
  d->CoordinateReferenceButtonGroup->addButton(d->LocalRadioButton, qMRMLTransformSliders::LOCAL);
  
  // Connect button group
  this->connect(d->CoordinateReferenceButtonGroup, SIGNAL(buttonPressed(int)), SLOT(onCoordinateReferenceButtonPressed(int)));

  // Connect identity button
  this->connect(d->IdentityPushButton, SIGNAL(clicked()), SLOT(identity()));

  // Connect revert button
  this->connect(d->InvertPushButton, SIGNAL(clicked()), SLOT(invert()));

  // Connect node selector with module itself
  this->connect(d->TransformNodeSelector, SIGNAL(currentNodeChanged(vtkMRMLNode*)), SLOT(onNodeSelected(vtkMRMLNode*)));

  // Connect minimum and maximum from the translation sliders to the matrix
  this->connect(d->TranslationSliders, SIGNAL(rangeChanged(double,double)), SLOT(onTranslationRangeChanged(double,double)));

  // Notify the matrix of the current translation min/max values
  this->onTranslationRangeChanged(d->TranslationSliders->minimum(), d->TranslationSliders->maximum());

  // Transform nodes connection
  this->connect(d->TransformToolButton, SIGNAL(clicked()), SLOT(transformSelectedNodes())); 
  this->connect(d->UntransformToolButton, SIGNAL(clicked()), SLOT(untransformSelectedNodes()));

  // Icons
  QIcon rightIcon = QApplication::style()->standardIcon(QStyle::SP_ArrowRight);
  d->TransformToolButton->setIcon(rightIcon);

  QIcon leftIcon = QApplication::style()->standardIcon(QStyle::SP_ArrowLeft);
  d->UntransformToolButton->setIcon(leftIcon);

  this->onNodeSelected(0);
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onFiducalNodeChanged(vtkMRMLNode* node)
{
	std::vector<vtkMRMLNode*> nodes;
	
	this->mrmlScene()->GetNodesByClass("vtkMRMLAnnotationHierarchyNode", nodes);
				
	if (nodes.size() > 1 && (d_ptr->inputFiducialsNodeSelector->currentNode() != 0))
	{		
		for (unsigned int i=0; i<nodes.size(); i++)
		{			
			vtkSmartPointer<vtkCollection> cnodes=vtkSmartPointer<vtkCollection>::New();
			vtkMRMLAnnotationHierarchyNode::SafeDownCast(nodes[i])->GetDirectChildren(cnodes);
			
			if (cnodes->GetNumberOfItems() > 1 && vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0)) != NULL 
				&& vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1)) != NULL)
			{
				fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0));
				fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1));
				
				if(!observersAreActive)
				{
					CallBack1->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial1Moved);
					CallBack2->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial2Moved);
					
					fnode0->AddObserver(fnode0->ControlPointModifiedEvent, CallBack1);
					fnode1->AddObserver(fnode1->ControlPointModifiedEvent, CallBack2);
					observersAreActive=true;
				}

				fnode0->GetFiducialCoordinates(startPoint); 
				fnode1->GetFiducialCoordinates(targetPoint);				
			}
		}
	}
	else
	{
		std::vector<vtkMRMLNode*> modelNodes;
		this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "ApplicatorPath")==0) // There is already a model:Applicator update its polydata
			{			
				//this->mrmlScene()->RemoveNode(this->mrmlScene()->GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
				//this->mrmlScene()->RemoveNode(modelNodes[j]);
				break;
			}
		}	
		this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);

		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "Applicator")==0) // There is already a model:Applicator update its polydata
			{			
				//this->mrmlScene()->RemoveNode(this->mrmlScene()->GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
				//his->mrmlScene()->RemoveNode(modelNodes[j]);
				break;
			}
		}
		modelNodes.clear();
	}
	nodes.clear();		
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::CreatePath()
{	
	std::vector<vtkMRMLNode*> nodes;	
	this->mrmlScene()->GetNodesByClass("vtkMRMLAnnotationHierarchyNode", nodes);
				
	if (nodes.size() < 2 || d_ptr->inputFiducialsNodeSelector->currentNode() == 0) // At least 2 nodes must exist(fiducials to create path)
	{
		nodes.clear(); 
		return;
	}

	for (unsigned int i=0; i<nodes.size(); i++)
	{			
		vtkSmartPointer<vtkCollection> cnodes=vtkSmartPointer<vtkCollection>::New();
		vtkMRMLAnnotationHierarchyNode::SafeDownCast(nodes[i])->GetDirectChildren(cnodes);
		
		if (cnodes->GetNumberOfItems() > 1 && vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0)) != NULL 
			&& vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1)) != NULL) // There are at least two fiducial child nodes
		{
			fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0));
			fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1));

			if(fnode0->GetFiducialCoordinates(startPoint) && fnode1->GetFiducialCoordinates(targetPoint))
			{
				CallBack1->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial1Moved);
				CallBack2->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial2Moved);

				fnode0->AddObserver(fnode0->ControlPointModifiedEvent, CallBack1);
				fnode1->AddObserver(fnode1->ControlPointModifiedEvent, CallBack2);

				std::vector<vtkMRMLNode*> modelNodes;
				this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
				
				// Define a line for the path
				vtkSmartPointer<vtkLineSource> line=vtkSmartPointer<vtkLineSource>::New();
				vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();						
				line->SetPoint1(startPoint[0], startPoint[1], startPoint[2]);
				line->SetPoint2(targetPoint[0], targetPoint[1], targetPoint[2]);									
				polydata=line->GetOutput();

				// Define a cylinder for the applicator
				//vtkSmartPointer<vtkCylinderSource> cylinder=vtkSmartPointer<vtkCylinderSource>::New(); // Create the applicator as a cylinder				
				//double x=targetPoint[0]-startPoint[0]; // Set the direction vector of the applicator
				//double y=targetPoint[1]-startPoint[1];
				//double z=targetPoint[2]-startPoint[2];
				//double h=50; // Height of cylinder
				//double r=2;  // Radius of cylinder
				//double angle=(180/3.1415926)*acos(y/sqrt(x*x+y*y+z*z)); // The angle between the applicator and the y-axis (default axis of creation)
				//cylinder->SetHeight(h); 
				//cylinder->SetRadius(r); 
				//cylinder->SetCenter(0.0, h/2., 0.0); 
				//cylinder->SetResolution(16);								
				
				// Define a transform for the cylinder
				vtkSmartPointer<vtkTransform> transform=vtkSmartPointer<vtkTransform>::New();
				//transform->Translate(startPoint[0], startPoint[1], startPoint[2]);
				//transform->RotateWXYZ(angle, z, 0, -1*x);				
				//vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();								
				//polydata=cylinder->GetOutput();

				vtkSmartPointer<vtkMRMLLinearTransformNode> transformNode=vtkSmartPointer<vtkMRMLLinearTransformNode>::New();
				transformNode->SetAndObserveMatrixTransformToParent( transform->GetMatrix() );
				transformNode->SetName("Transform Path");


				bool found=false;
									
				for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
				{
					if(strcmp(modelNodes[j]->GetName(), "ApplicatorPath")==0) // There is already a model:Applicator update its polydata
					{	
						this->mrmlScene()->RemoveNode(this->mrmlScene()->
							GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
						this->mrmlScene()->AddNode(transformNode);
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObserveTransformNodeID(transformNode->GetID());
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(1);
											
						pathVisible=true;
						found=true;
						break;
					}
				}				

				if(!found)
				{
					this->mrmlScene()->AddNode(transformNode);
					vtkSmartPointer<vtkMRMLModelNode> modelNode=vtkSmartPointer<vtkMRMLModelNode>::New();							
					modelNode->SetAndObserveTransformNodeID(transformNode->GetID());
					modelNode->SetScene(this->mrmlScene());
					modelNode->SetName("ApplicatorPath");
					modelNode->SetAndObservePolyData(polydata);					

					vtkSmartPointer<vtkMRMLModelDisplayNode> modelDisplayNode=vtkSmartPointer<vtkMRMLModelDisplayNode>::New();							
					modelDisplayNode->SetColor(1,0,0); // Red
					modelDisplayNode->SetVisibility(1);
					pathVisible=true;

					this->mrmlScene()->AddNode(modelDisplayNode);
					modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

					this->mrmlScene()->AddNode(modelNode);	
				}
				modelNodes.clear();
			}
		}
	}	
	nodes.clear();	
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::HidePath()
{
	std::vector<vtkMRMLNode*> modelNodes;
	this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
	for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
	{
		if(strcmp(modelNodes[j]->GetName(), "ApplicatorPath")==0) // There is already a model:Applicator update its polydata
		{			
			this->mrmlScene()->RemoveNode(this->mrmlScene()->GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
			this->mrmlScene()->RemoveNode(modelNodes[j]);
			break;
		}
	}	
	modelNodes.clear();
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::LoadApplicator()
{	
	/*d_ptr->btnPath->setText(QString::number(count));
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("C:\\Users\\eyeniaras\\Downloads\\daviddata\\laserApplicator.stl");*/
	
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

	//vtkMRMLLinearTransformNode* transformNode=vtkMRMLLinearTransformNode::New();
	//this->mrmlScene()->AddNode(transformNode);

	//vtkSmartPointer<vtkMatrix4x4> mat = vtkSmartPointer<vtkMatrix4x4>::New();
 //   mat->Identity();
	//double x=targetPoint[0]-startPoint[0];
	//double y=targetPoint[1]-startPoint[1];
	//double z=targetPoint[2]-startPoint[2];
	//double h=sqrt(x*x+y*y+z*z);
	//double angle=acos(y/h);
	//angle=angle*(180/3.1415926); // Convert from radians to degrees

	////double h1=sqrt(x*x+y*y);
	////double h2=sqrt(z*z+y*y);
	////double cosT=y/h2, sinT=z/h2, cosA=y/h1, sinA=x/h1;

	///*mat->SetElement(0, 0, cosA);      mat->SetElement(0, 1, -1*sinA);   mat->SetElement(0, 2, 0);
	//mat->SetElement(1, 0, cosT*sinA); mat->SetElement(1, 1, cosT*cosA); mat->SetElement(1, 2, -1*sinT);
	//mat->SetElement(2, 0, sinT*sinA); mat->SetElement(2, 1, sinT*cosA); mat->SetElement(2, 2, cosT);		
	//mat->SetElement(0,3,startPoint[0]);
	//mat->SetElement(1,3,startPoint[1]);
	//mat->SetElement(2,3,startPoint[2]);
	//transformNode->SetAndObserveMatrixTransformToParent( mat.GetPointer() );
	//*/

	//vtkSmartPointer<vtkTransform> transform=vtkSmartPointer<vtkTransform>::New();
	//
	//transform->Translate(startPoint[0] + x/2, startPoint[1] + y/2, startPoint[2] + z/2);
	//transform->RotateWXYZ(angle, z, 0, -1*x);

	////transform->RotateX(atan(z/y)*(180/3.14));
	//d_ptr->btnBurn->setText(QString::number(angle));
	////transformNode->ApplyTransform(transform);

	////transformNode->ApplyTransformMatrix(transform->GetMatrix());
	//transformNode->SetAndObserveMatrixTransformToParent( transform->GetMatrix() );
	//
	//vtkMRMLModelNode* modelNode=vtkMRMLModelNode::New();
	//modelNode->SetAndObserveTransformNodeID(transformNode->GetID());
	//modelNode->SetScene(this->mrmlScene());
	//modelNode->SetName("Applicator");
	//modelNode->SetAndObservePolyData(reader->GetOutput());	
 //   reader->Update();

	//vtkMRMLModelDisplayNode* modelDisplayNode=vtkMRMLModelDisplayNode::New();
	//modelDisplayNode->SetColor(1,0,0); // green
	//modelDisplayNode->SetScene(this->mrmlScene());
	//this->mrmlScene()->AddNode(modelDisplayNode);
	//modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

	//modelDisplayNode->SetPolyData(modelNode->GetPolyData());
	//this->mrmlScene()->AddNode(modelNode);

	//transformNode->Delete();
	//modelNode->Delete();
	//modelDisplayNode->Delete();	    

	std::vector<vtkMRMLNode*> nodes;	
	this->mrmlScene()->GetNodesByClass("vtkMRMLAnnotationHierarchyNode", nodes);
				
	if (nodes.size() < 2 || d_ptr->inputFiducialsNodeSelector->currentNode() == 0) // At least 2 nodes must exist(fiducials to create path)
	{
		nodes.clear(); 
		return;
	}

	for (unsigned int i=0; i<nodes.size(); i++)
	{			
		vtkSmartPointer<vtkCollection> cnodes=vtkSmartPointer<vtkCollection>::New();
		vtkMRMLAnnotationHierarchyNode::SafeDownCast(nodes[i])->GetDirectChildren(cnodes);
		
		if (cnodes->GetNumberOfItems() > 1 && vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0)) != NULL 
			&& vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1)) != NULL) // There are at least two fiducial child nodes
		{
			fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0));
			fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1));

			if(fnode0->GetFiducialCoordinates(startPoint) && fnode1->GetFiducialCoordinates(targetPoint))
			{
				CallBack1->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial1Moved);
				CallBack2->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial2Moved);

				fnode0->AddObserver(fnode0->ControlPointModifiedEvent, CallBack1);
				fnode1->AddObserver(fnode1->ControlPointModifiedEvent, CallBack2);

				std::vector<vtkMRMLNode*> modelNodes;
				this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
				
				// Define a cylinder for the applicator
				vtkSmartPointer<vtkCylinderSource> cylinder=vtkSmartPointer<vtkCylinderSource>::New(); // Create the applicator as a cylinder				
				double x=targetPoint[0]-startPoint[0]; // Set the direction vector of the applicator
				double y=targetPoint[1]-startPoint[1];
				double z=targetPoint[2]-startPoint[2];
				double h=50; // Height of cylinder
				double r=2;  // Radius of cylinder
				double angle=(180/3.1415926)*acos(y/sqrt(x*x+y*y+z*z)); // The angle between the applicator and the y-axis (default axis of creation)
				cylinder->SetHeight(h); 
				cylinder->SetRadius(r); 
				cylinder->SetCenter(0.0, h/2., 0.0); 
				cylinder->SetResolution(16);								
				
				// Define a transform for the cylinder
				vtkSmartPointer<vtkTransform> transform=vtkSmartPointer<vtkTransform>::New();
				transform->Translate(startPoint[0], startPoint[1], startPoint[2]);
				transform->RotateWXYZ(angle, z, 0, -1*x);				
				vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();								
				polydata=cylinder->GetOutput();

				vtkSmartPointer<vtkMRMLLinearTransformNode> transformNode=vtkSmartPointer<vtkMRMLLinearTransformNode>::New();
				transformNode->SetAndObserveMatrixTransformToParent( transform->GetMatrix() );
				transformNode->SetName("Transform Applicator");

				bool found=false;
									
				for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
				{
					if(strcmp(modelNodes[j]->GetName(), "Applicator")==0) // There is already a model:Applicator update its polydata
					{	
						this->mrmlScene()->RemoveNode(this->mrmlScene()->
							GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
						this->mrmlScene()->AddNode(transformNode);
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObserveTransformNodeID(transformNode->GetID());
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(1);
											
						applicatorVisible=true;
						found=true;
						break;
					}
				}				

				if(!found)
				{
					this->mrmlScene()->AddNode(transformNode);
					vtkSmartPointer<vtkMRMLModelNode> modelNode=vtkSmartPointer<vtkMRMLModelNode>::New();							
					modelNode->SetAndObserveTransformNodeID(transformNode->GetID());
					modelNode->SetScene(this->mrmlScene());
					modelNode->SetName("Applicator");
					modelNode->SetAndObservePolyData(polydata);					

					vtkSmartPointer<vtkMRMLModelDisplayNode> modelDisplayNode=vtkSmartPointer<vtkMRMLModelDisplayNode>::New();							
					modelDisplayNode->SetColor(0,1,0); // Red
					modelDisplayNode->SetVisibility(1);
					applicatorVisible=true;

					this->mrmlScene()->AddNode(modelDisplayNode);
					modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

					this->mrmlScene()->AddNode(modelNode);	
				}
				modelNodes.clear();
			}
		}
	}	
	nodes.clear();	
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::CreateSphereAtTarget()
{	
	std::vector<vtkMRMLNode*> nodes;	
	this->mrmlScene()->GetNodesByClass("vtkMRMLAnnotationHierarchyNode", nodes);
				
	if (nodes.size() < 2 || d_ptr->inputFiducialsNodeSelector->currentNode() == 0) // At least 2 nodes must exist(fiducials to create path)
	{
		nodes.clear(); 
		return;
	}

	for (unsigned int i=0; i<nodes.size(); i++)
	{			
		vtkSmartPointer<vtkCollection> cnodes=vtkSmartPointer<vtkCollection>::New();
		vtkMRMLAnnotationHierarchyNode::SafeDownCast(nodes[i])->GetDirectChildren(cnodes);
		
		if (cnodes->GetNumberOfItems() > 1 && vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0)) != NULL 
			&& vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1)) != NULL) // There are at least two fiducial child nodes
		{
			fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0));
			fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1));

			if(fnode0->GetFiducialCoordinates(startPoint) && fnode1->GetFiducialCoordinates(targetPoint))
			{
				CallBack1->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial1Moved);
				CallBack2->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial2Moved);

				fnode0->AddObserver(fnode0->ControlPointModifiedEvent, CallBack1);
				fnode1->AddObserver(fnode1->ControlPointModifiedEvent, CallBack2);

				std::vector<vtkMRMLNode*> modelNodes;
				this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
				
				// Define a cylinder for the applicator
				vtkSmartPointer<vtkSphereSource> sphere=vtkSmartPointer<vtkSphereSource>::New(); 	
				vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();
				sphere->SetCenter(targetPoint[0], targetPoint[1], targetPoint[2]);
				sphere->SetRadius(5); 	
				
				sphere->SetThetaResolution(32);
				sphere->SetPhiResolution(32);
				sphere->SetStartTheta(0.0);
				sphere->SetEndTheta(360.0);
				sphere->SetStartPhi(0.0);
				sphere->SetEndPhi(360.0);
				sphere->SetLatLongTessellation(1);
				polydata=sphere->GetOutput();

				sleep(500);

				
					
				
				// Define a transform for the cylinder
				vtkSmartPointer<vtkTransform> transform=vtkSmartPointer<vtkTransform>::New();
				
				vtkSmartPointer<vtkMRMLLinearTransformNode> transformNode=vtkSmartPointer<vtkMRMLLinearTransformNode>::New();
				transformNode->SetAndObserveMatrixTransformToParent( transform->GetMatrix() );
				transformNode->SetName("Transform Treatment");


				bool found=false;
									
				for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
				{
					if(strcmp(modelNodes[j]->GetName(), "Treatment")==0) // There is already a model:Applicator update its polydata
					{	
						this->mrmlScene()->RemoveNode(this->mrmlScene()->
							GetNodeByID(vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetTransformNodeID()));
						this->mrmlScene()->AddNode(transformNode);
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObserveTransformNodeID(transformNode->GetID());
						vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(1);
											
						//pathVisible=true;
						found=true;
						break;
					}
				}				

				if(!found)
				{
					this->mrmlScene()->AddNode(transformNode);
					vtkSmartPointer<vtkMRMLModelNode> modelNode=vtkSmartPointer<vtkMRMLModelNode>::New();							
					modelNode->SetAndObserveTransformNodeID(transformNode->GetID());
					modelNode->SetScene(this->mrmlScene());
					modelNode->SetName("Treatment");
					modelNode->SetAndObservePolyData(polydata);					

					vtkSmartPointer<vtkMRMLModelDisplayNode> modelDisplayNode=vtkSmartPointer<vtkMRMLModelDisplayNode>::New();							
					modelDisplayNode->SetColor(0,0,1); // Red
					modelDisplayNode->SetVisibility(1);
					//pathVisible=true;

					this->mrmlScene()->AddNode(modelDisplayNode);
					modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

					this->mrmlScene()->AddNode(modelNode);	
				}
				modelNodes.clear();
			}
		}
	}	
	nodes.clear();	
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::OnFiducial1Moved(vtkObject* vtk_obj, unsigned long event, void* client_data, void* call_data)
{
	vtkMRMLAnnotationFiducialNode* fnode = reinterpret_cast<vtkMRMLAnnotationFiducialNode*>(vtk_obj);
	qSlicerLITTPlanV2ModuleWidget* thisClass = reinterpret_cast<qSlicerLITTPlanV2ModuleWidget*>(client_data);
	
	//if(thisClass->applicatorVisible || thisClass->pathVisible)
	{
		fnode->GetFiducialCoordinates(thisClass->startPoint);

		/*vtkSmartPointer<vtkLineSource> line=vtkSmartPointer<vtkLineSource>::New();
		vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();						
		line->SetPoint1(thisClass->startPoint[0], thisClass->startPoint[1], thisClass->startPoint[2]);
		line->SetPoint2(thisClass->targetPoint[0], thisClass->targetPoint[1], thisClass->targetPoint[2]);									
		polydata=line->GetOutput();*/

		std::vector<vtkMRMLNode*> modelNodes;	
		thisClass->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
		
		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "ApplicatorPath")==0) // There is already a model:ApplicatorPath update its polydata
			{	
				//if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
				//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
				//thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
				vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
				//thisClass->applicatorVisible=false;
				break;
			}
		}

		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "Applicator")==0) // There is already a model:ApplicatorPath update its polydata
			{	
				//if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
				//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
				//thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
				vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
				//thisClass->pathVisible=false;
				break;
			}
		}
		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "Treatment")==0) // There is already a model:ApplicatorPath update its polydata
			{	
				//if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
				//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
				//thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
				vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
				//thisClass->pathVisible=false;
				break;
			}
		}
		modelNodes.clear();
	}

	/*for(unsigned int c=0; c<3; c++)
	{
		thisClass->prevStartPoint[c]=thisClass->startPoint[c];
	}*/

	//thisClass->CreatePath();
	
	//vtkMRMLAnnotationFiducialNode* fnode0 = reinterpret_cast<vtkMRMLAnnotationFiducialNode*>(vtk_obj);

	//if(fnode0->GetFiducialCoordinates(qSlicerLITTPlanV2ModuleWidget::startPoint))
	//{
	//}

	//vtkMRMLAnnotationFiducialNode *fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(vtk_obj);
	//vtkMRMLAnnotationFiducialNode *fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(vtk_obj);
	
	
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::OnFiducial2Moved(vtkObject* vtk_obj, unsigned long event, void* client_data, void* call_data)
{
	vtkMRMLAnnotationFiducialNode* fnode = reinterpret_cast<vtkMRMLAnnotationFiducialNode*>(vtk_obj);
	qSlicerLITTPlanV2ModuleWidget* thisClass = reinterpret_cast<qSlicerLITTPlanV2ModuleWidget*>(client_data);
	
	//if(thisClass->applicatorVisible || thisClass->pathVisible)
	{
		fnode->GetFiducialCoordinates(thisClass->targetPoint);

		//vtkSmartPointer<vtkLineSource> line=vtkSmartPointer<vtkLineSource>::New();
		//vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();						
		//line->SetPoint1(thisClass->startPoint[0], thisClass->startPoint[1], thisClass->startPoint[2]);
		//line->SetPoint2(thisClass->targetPoint[0], thisClass->targetPoint[1], thisClass->targetPoint[2]);									
		//polydata=line->GetOutput();

		std::vector<vtkMRMLNode*> modelNodes;	
		thisClass->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
		
		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "ApplicatorPath")==0) // There is already a model:ApplicatorPath update its polydata
			{	
				//if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
				//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
				//thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
				vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
				//thisClass->applicatorVisible=false;
				break;
			}
		}

		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "Applicator")==0) // There is already a model:ApplicatorPath update its polydata
			{	
				//if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
				//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
				//thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
				vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
				//thisClass->pathVisible=false;
				break;
			}
		}

		for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
		{
			if(strcmp(modelNodes[j]->GetName(), "Treatment")==0) // There is already a model:ApplicatorPath update its polydata
			{	
				//if(fabs(thisClass->startPoint[0]-thisClass->prevStartPoint[0])>20)
				//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
				//thisClass->d_ptr->btnPath->setText(QString::number(thisClass->startPoint[0]));
				vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(0);
				//thisClass->pathVisible=false;
				break;
			}
		}
		modelNodes.clear();
	}

	/*for(unsigned int c=0; c<3; c++)
	{
		thisClass->prevStartPoint[c]=thisClass->startPoint[c];
	}*/

	//thisClass->CreatePath();
	
	//vtkMRMLAnnotationFiducialNode* fnode0 = reinterpret_cast<vtkMRMLAnnotationFiducialNode*>(vtk_obj);

	//if(fnode0->GetFiducialCoordinates(qSlicerLITTPlanV2ModuleWidget::startPoint))
	//{
	//}

	//vtkMRMLAnnotationFiducialNode *fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(vtk_obj);
	//vtkMRMLAnnotationFiducialNode *fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(vtk_obj);
	
	
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onTissueTypeChanged()
{
	LoadParametersFromGUI(); // First get the arameter values from screen (the values before the new tissue selection)
	LoadParametersToGUI();	// Now load the selected tissue`s parameters to the screen
	tissueIndex=d_ptr->TissueSelector->currentIndex();
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::SetParametersToDefault()
{
	this->wmtParameters[0]=1.1;
	this->wmtParameters[1]=1.2;
	this->wmtParameters[2]=1.3;
	this->wmtParameters[3]=1.4;
	this->wmtParameters[4]=1.5;

	this->gmtParameters[0]=2.1;
	this->gmtParameters[1]=2.2;
	this->gmtParameters[2]=2.3;
	this->gmtParameters[3]=2.3;
	this->gmtParameters[4]=2.5;

	this->csfParameters[0]=3.1;
	this->csfParameters[1]=3.1;
	this->csfParameters[2]=3.1;
	this->csfParameters[3]=3.1;
	this->csfParameters[4]=3.1;

	this->tumParameters[0]=4.0;
	this->tumParameters[1]=4.0;
	this->tumParameters[2]=4.0;
	this->tumParameters[3]=4.0;
	this->tumParameters[4]=4.0;	

	power=5.00;
	thermalDose=5.00;
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::LoadParametersToGUI()
{
	switch(d_ptr->TissueSelector->currentIndex())
	{
	case 0:
		d_ptr->ThermalConductivity->setText("0.00");
		d_ptr->TissuePerfusion->setText("0.00");
		d_ptr->OpticalAbsorption->setText("0.00");
		d_ptr->OpticalScattering->setText("0.00");
		d_ptr->OpticalAnisotrophy->setText("0.00"); 
		break;
	case 1:
		d_ptr->ThermalConductivity->setText(QString::number(wmtParameters[0]));
		d_ptr->TissuePerfusion->setText(QString::number(wmtParameters[1]));
		d_ptr->OpticalAbsorption->setText(QString::number(wmtParameters[2]));
		d_ptr->OpticalScattering->setText(QString::number(wmtParameters[3]));
		d_ptr->OpticalAnisotrophy->setText(QString::number(wmtParameters[4])); 
		break;
	case 2:
		d_ptr->ThermalConductivity->setText(QString::number(gmtParameters[0]));
		d_ptr->TissuePerfusion->setText(QString::number(gmtParameters[1]));
		d_ptr->OpticalAbsorption->setText(QString::number(gmtParameters[2]));
		d_ptr->OpticalScattering->setText(QString::number(gmtParameters[3]));
		d_ptr->OpticalAnisotrophy->setText(QString::number(gmtParameters[4])); 
		break;
	case 3:
		d_ptr->ThermalConductivity->setText(QString::number(csfParameters[0]));
		d_ptr->TissuePerfusion->setText(QString::number(csfParameters[1]));
		d_ptr->OpticalAbsorption->setText(QString::number(csfParameters[2]));
		d_ptr->OpticalScattering->setText(QString::number(csfParameters[3]));
		d_ptr->OpticalAnisotrophy->setText(QString::number(csfParameters[4])); 
		break;
	case 4:
		d_ptr->ThermalConductivity->setText(QString::number(tumParameters[0]));
		d_ptr->TissuePerfusion->setText(QString::number(tumParameters[1]));
		d_ptr->OpticalAbsorption->setText(QString::number(tumParameters[2]));
		d_ptr->OpticalScattering->setText(QString::number(tumParameters[3]));
		d_ptr->OpticalAnisotrophy->setText(QString::number(tumParameters[4]));  
		break;
	default:
		break;
	}
	
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::LoadParametersFromGUI()
{
	switch(tissueIndex)
	{
	case 0:
		break;
	case 1:
		this->wmtParameters[0]=d_ptr->ThermalConductivity->text().toDouble();
		this->wmtParameters[1]=d_ptr->TissuePerfusion->text().toDouble();
		this->wmtParameters[2]=d_ptr->OpticalAbsorption->text().toDouble();
		this->wmtParameters[3]=d_ptr->OpticalScattering->text().toDouble();
		this->wmtParameters[4]=d_ptr->OpticalAnisotrophy->text().toDouble();
		break;
	case 2:
		this->gmtParameters[0]=d_ptr->ThermalConductivity->text().toDouble();
		this->gmtParameters[1]=d_ptr->TissuePerfusion->text().toDouble();
		this->gmtParameters[2]=d_ptr->OpticalAbsorption->text().toDouble();
		this->gmtParameters[3]=d_ptr->OpticalScattering->text().toDouble();
		this->gmtParameters[4]=d_ptr->OpticalAnisotrophy->text().toDouble();
		break;
	case 3:
		this->csfParameters[0]=d_ptr->ThermalConductivity->text().toDouble();
		this->csfParameters[1]=d_ptr->TissuePerfusion->text().toDouble();
		this->csfParameters[2]=d_ptr->OpticalAbsorption->text().toDouble();
		this->csfParameters[3]=d_ptr->OpticalScattering->text().toDouble();
		this->csfParameters[4]=d_ptr->OpticalAnisotrophy->text().toDouble();
		break;
	case 4:
		this->tumParameters[0]=d_ptr->ThermalConductivity->text().toDouble();
		this->tumParameters[1]=d_ptr->TissuePerfusion->text().toDouble();
		this->tumParameters[2]=d_ptr->OpticalAbsorption->text().toDouble();
		this->tumParameters[3]=d_ptr->OpticalScattering->text().toDouble();
		this->tumParameters[4]=d_ptr->OpticalAnisotrophy->text().toDouble();
		break;
	default:
		break;
	}
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::onBtnDefaultClicked()
{
	SetParametersToDefault();
	LoadParametersToGUI();
	d_ptr->ThermalDose->setText(QString::number(thermalDose));
	d_ptr->Power->setText(QString::number(power));
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
void qSlicerLITTPlanV2ModuleWidget::onBtnBurnClicked() // Export the created points here!! FEM will be read too..
{
	// Update all the parameters from GUI
	LoadParametersFromGUI(); // Tissue parameters
	this->thermalDose=d_ptr->ThermalDose->text().toDouble(); // Global parameter: Thermal dose
	this->power=d_ptr->Power->text().toDouble(); // Global parameter: power
		
	// Write all the parameter values and the two guiding points to .ini file
	ofstream outputFile;
	outputFile.open("LITTData.ini");
	
	outputFile<<"[EntrancePoint]"<<endl;
	outputFile<<"x="<<startPoint[0]<<endl;
	outputFile<<"y="<<startPoint[1]<<endl;
	outputFile<<"z="<<startPoint[2]<<endl;
	
	outputFile<<"[TargetPoint]"<<endl;
	outputFile<<"x="<<targetPoint[0]<<endl;
	outputFile<<"y="<<targetPoint[1]<<endl;
	outputFile<<"z="<<targetPoint[2]<<endl;

	outputFile<<"[WhiteMatter]"<<endl;
	outputFile<<"ThermalConductivity="<<wmtParameters[0]<<endl;
	outputFile<<"TissuePerfusion="<<wmtParameters[1]<<endl;
	outputFile<<"OpticalAbsorbtion="<<wmtParameters[2]<<endl;
	outputFile<<"OpticalScattering="<<wmtParameters[3]<<endl;
	outputFile<<"OpticalAnisotrophy="<<wmtParameters[4]<<endl;

	outputFile<<"[GreyMatter]"<<endl;
	outputFile<<"ThermalConductivity="<<gmtParameters[0]<<endl;
	outputFile<<"TissuePerfusion="<<gmtParameters[1]<<endl;
	outputFile<<"OpticalAbsorbtion="<<gmtParameters[2]<<endl;
	outputFile<<"OpticalScattering="<<gmtParameters[3]<<endl;
	outputFile<<"OpticalAnisotrophy="<<gmtParameters[4]<<endl;

	outputFile<<"[CSF]"<<endl;
	outputFile<<"ThermalConductivity="<<csfParameters[0]<<endl;
	outputFile<<"TissuePerfusion="<<csfParameters[1]<<endl;
	outputFile<<"OpticalAbsorbtion="<<csfParameters[2]<<endl;
	outputFile<<"OpticalScattering="<<csfParameters[3]<<endl;
	outputFile<<"OpticalAnisotrophy="<<csfParameters[4]<<endl;

	outputFile<<"[Tumor]"<<endl;
	outputFile<<"ThermalConductivity="<<tumParameters[0]<<endl;
	outputFile<<"TissuePerfusion="<<tumParameters[1]<<endl;
	outputFile<<"OpticalAbsorbtion="<<tumParameters[2]<<endl;
	outputFile<<"OpticalScattering="<<tumParameters[3]<<endl;
	outputFile<<"OpticalAnisotrophy="<<tumParameters[4]<<endl;

	outputFile<<"[Global]"<<endl;
	outputFile<<"ThermalDose="<<thermalDose<<endl;
	outputFile<<"Power="<<power<<endl;

	outputFile.close();

	CreateSphereAtTarget();	
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
void qSlicerLITTPlanV2ModuleWidget::extractMinMaxTranslationValue(vtkMatrix4x4 * mat, double& min, double& max)
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
void qSlicerLITTPlanV2ModuleWidget::onTranslationRangeChanged(double newMin, double newMax)
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

//-----------------------------------------------------------------------------
bool qSlicerLITTPlanV2ModuleWidget::EqualD(double a, double b)
{
	double dTolerance=0.000001;

	if(fabs(a-b)<dTolerance) 
		return true;
	else 
		return false;
}
//-----------------------------------------------------------------------------
bool qSlicerLITTPlanV2ModuleWidget::EqualP(double p1[], double p2[])
{
	if( EqualD(p1[0], p2[0]) && EqualD(p1[1], p2[1])&& EqualD(p1[2], p2[2]))
		return true;
	else 
		return false;
}

//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::CreatePathOld()
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
	std::vector<vtkMRMLNode*> modelNodes;

	this->mrmlScene()->GetNodesByClass("vtkMRMLAnnotationHierarchyNode", nodes);
	this->mrmlScene()->GetNodesByClass("vtkMRMLModelNode", modelNodes);
			
	if (nodes.size() > 1 && (d_ptr->inputFiducialsNodeSelector->currentNode() != 0))
	{		
		for (unsigned int i=0; i<nodes.size(); i++)
		{			
			vtkSmartPointer<vtkCollection> cnodes=vtkSmartPointer<vtkCollection>::New();
			vtkMRMLAnnotationHierarchyNode::SafeDownCast(nodes[i])->GetDirectChildren(cnodes);
			
			if (cnodes->GetNumberOfItems() > 1 && vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0)) != NULL 
				&& vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1)) != NULL)
			{
				fnode0 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(0));
				fnode1 = vtkMRMLAnnotationFiducialNode::SafeDownCast(cnodes->GetItemAsObject(1));
				
				if(!observersAreActive)
				{
				
					CallBack1->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial1Moved);
					CallBack2->SetCallback(qSlicerLITTPlanV2ModuleWidget::OnFiducial2Moved);
				
					fnode0->AddObserver(fnode0->ControlPointModifiedEvent, CallBack1);
					fnode1->AddObserver(fnode1->ControlPointModifiedEvent, CallBack2);
					observersAreActive=true;
				}
				
				if(fnode0->GetFiducialCoordinates(startPoint) && fnode1->GetFiducialCoordinates(targetPoint))
				{
					fiducialsAreValid=true;
					//if( (!EqualP(startPoint, prevStartPoint)) || (!EqualP(targetPoint, prevTargetPoint))) // If any of the points has been changed (=different than the previous value)
					{						
						vtkSmartPointer<vtkLineSource> line=vtkSmartPointer<vtkLineSource>::New();
						vtkSmartPointer<vtkPolyData> polydata=vtkSmartPointer<vtkPolyData>::New();						
						line->SetPoint1(startPoint[0], startPoint[1], startPoint[2]);
						line->SetPoint2(targetPoint[0], targetPoint[1], targetPoint[2]);									
						polydata=line->GetOutput();

						bool found=false;
						
						for(unsigned int j=0; j<modelNodes.size(); j++) // Check all the models
						{
							if(strcmp(modelNodes[j]->GetName(), "ApplicatorPath")==0) // There is already a model:ApplicatorPath update its polydata
							{
								vtkMRMLModelNode* pathModel=vtkMRMLModelNode::SafeDownCast(modelNodes[j]);
								pathModel->SetAndObservePolyData(polydata);
								pathModel->GetModelDisplayNode()->SetVisibility(1);
								pathVisible=true;
								//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->SetAndObservePolyData(polydata);
								//vtkMRMLModelNode::SafeDownCast(modelNodes[j])->GetModelDisplayNode()->SetVisibility(1);
								found=true;
								break;
							}
						}
						
						if(!found) // create a new one
						{
							vtkSmartPointer<vtkMRMLModelNode> modelNode=vtkSmartPointer<vtkMRMLModelNode>::New();							
							modelNode->SetScene(this->mrmlScene());
							modelNode->SetName("ApplicatorPath");							
							modelNode->SetAndObservePolyData(polydata);	
												
							vtkSmartPointer<vtkMRMLModelDisplayNode> modelDisplayNode=vtkSmartPointer<vtkMRMLModelDisplayNode>::New();							
							modelDisplayNode->SetColor(1,0,0); // Red
							pathVisible=true;
							//modelDisplayNode->SetVisibility(1);

							this->mrmlScene()->AddNode(modelDisplayNode);
							modelNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());

							this->mrmlScene()->AddNode(modelNode);													
						}				
						
						/*for(int c=0; c<3; c++)
						{
							prevStartPoint[c]=startPoint[c];
							prevTargetPoint[c]=targetPoint[c];
						}*/
					}
				}															
				break;
			}
		}
	}
	nodes.clear();	
	modelNodes.clear();
	//d_ptr->btnPath->setText(QString::number(pathVisible));
}
//-----------------------------------------------------------------------------
void qSlicerLITTPlanV2ModuleWidget::sleep(unsigned int mseconds)
{
    clock_t goal = mseconds + clock();
    while (goal > clock());
}