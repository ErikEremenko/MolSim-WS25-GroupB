/*
 * VTKWriter.cpp
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */
#ifdef ENABLE_VTK_OUTPUT

#include "io/VTKWriter.h"

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <filesystem>

#include <iomanip>
#include <sstream>

namespace outputWriter {

void VTKWriter::plotParticles(const ParticleContainer& particles, const std::string& filename, const int iteration) {
  // create separate output directory
  const std::string output_directory = "output";
  try {
    std::filesystem::create_directories(output_directory);
  } catch (const std::filesystem::filesystem_error& err) {
    std::cerr << "Error wile creating directory " << output_directory << ": " << err.what() << std::endl;
    return;
  }
  // Initialize points
  auto points = vtkSmartPointer<vtkPoints>::New();

  // Create and configure data arrays
  vtkNew<vtkFloatArray> massArray;
  massArray->SetName("mass");
  massArray->SetNumberOfComponents(1);

  vtkNew<vtkFloatArray> velocityArray;
  velocityArray->SetName("velocity");
  velocityArray->SetNumberOfComponents(3);

  vtkNew<vtkFloatArray> forceArray;
  forceArray->SetName("force");
  forceArray->SetNumberOfComponents(3);

  vtkNew<vtkIntArray> typeArray;
  typeArray->SetName("type");
  typeArray->SetNumberOfComponents(1);

  for (auto& p : particles) {
    points->InsertNextPoint(p.getX().data());
    massArray->InsertNextValue(static_cast<float>(p.getM()));
    velocityArray->InsertNextTuple(p.getV().data());
    forceArray->InsertNextTuple(p.getF().data());
    typeArray->InsertNextValue(p.getType());
  }

  // Set up the grid
  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->SetPoints(points);

  // Add arrays to the grid
  grid->GetPointData()->AddArray(massArray);
  grid->GetPointData()->AddArray(velocityArray);
  grid->GetPointData()->AddArray(forceArray);
  grid->GetPointData()->AddArray(typeArray);

  // Create filename with iteration number
  std::stringstream strstr;
  strstr << output_directory << "/" << filename << "_" << std::setfill('0') << std::setw(4) << iteration << ".vtu";

  // Create writer and set data
  vtkNew<vtkXMLUnstructuredGridWriter> writer;
  writer->SetFileName(strstr.str().c_str());
  writer->SetInputData(grid);
  writer->SetDataModeToAscii();
  writer->SetDataModeToBinary();

  // Write the file
  writer->Write();
  // TODO output to log
  /*std::cout << "Iteration " << iteration << " written." << std::endl;*/
}
}  // namespace outputWriter
#endif
