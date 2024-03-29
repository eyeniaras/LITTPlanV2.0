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

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

#ifndef __qSlicerLITTPlanV2Module_h
#define __qSlicerLITTPlanV2Module_h

// SlicerQt includes
#include "qSlicerLoadableModule.h"

// LITTPlanV2 includes
#include "qSlicerLITTPlanV2ModuleExport.h"

class vtkMatrix4x4;
class vtkMRMLNode;
class qSlicerLITTPlanV2ModulePrivate;

class Q_SLICER_QTMODULES_LITTPLANV2_EXPORT qSlicerLITTPlanV2Module
  : public qSlicerLoadableModule
{
  Q_OBJECT
  Q_INTERFACES(qSlicerLoadableModule);
public:

  typedef qSlicerLoadableModule Superclass;
  qSlicerLITTPlanV2Module(QObject *parent=0);
  virtual ~qSlicerLITTPlanV2Module();

  /// Icon of the transform module
  virtual QIcon icon()const;

  virtual QStringList categories()const;

  /// Display name for the module
  qSlicerGetTitleMacro("LITTPlanV2");

  /// Help text of the module
  virtual QString helpText()const;

  /// Acknowledgement of the module
  virtual QString acknowledgementText()const;

  /// Contributors of the module
  virtual QStringList contributors()const;

protected:
  /// Reimplemented to initialize the transforms IO
  virtual void setup();

  /// Create and return the widget representation associated to this module
  virtual qSlicerAbstractModuleRepresentation * createWidgetRepresentation();

  /// Create and return the logic associated to this module
  virtual vtkMRMLAbstractLogic* createLogic();

  QScopedPointer<qSlicerLITTPlanV2ModulePrivate> d_ptr;
private:
  Q_DECLARE_PRIVATE(qSlicerLITTPlanV2Module);
  Q_DISABLE_COPY(qSlicerLITTPlanV2Module);
};

#endif
