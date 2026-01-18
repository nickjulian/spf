import h5py
import vtk

def convert_phi_to_vti(file_path, step, output_file):
    with h5py.File(file_path, "r") as f:
        step_name = "Step{}".format(step)
        if step_name in f:
            dataset = f[step_name + "/phi"]

            nx, ny, nz = dataset.shape

            image_data = vtk.vtkImageData()
            image_data.SetDimensions(nx + 1, ny + 1, nz + 1)

            data_array = vtk.vtkFloatArray()
            data_array.SetNumberOfComponents(1)
            data_array.SetName("C")

            for i in range(nz):
                for j in range(ny):
                    for k in range(nx):
                        value = dataset[k, j, i]
                        data_array.InsertNextValue(value)

            image_data.GetCellData().SetScalars(data_array)
            image_data.SetSpacing(1.0, 1.0, 1.0)
            image_data.SetOrigin(0.0, 0.0, 0.0)

            writer = vtk.vtkXMLImageDataWriter()
            writer.SetCompressorTypeToNone()
            writer.SetDataModeToAscii()
            writer.SetFileName(output_file)
            writer.SetInputData(image_data)
            writer.Write()

def convert_T_to_vti(file_path, step, output_file):
    with h5py.File(file_path, "r") as f:
        step_name = "Step{}".format(step)
        if step_name in f:
            datasetT = f[step_name + "/T"]

            nx, ny, nz = datasetT.shape

            image_data = vtk.vtkImageData()
            image_data.SetDimensions(nx + 1, ny + 1, nz + 1)

            data_array = vtk.vtkFloatArray()
            data_array.SetNumberOfComponents(1)
            data_array.SetName("T")

            for i in range(nz):
                for j in range(ny):
                    for k in range(nx):
                        valueT = datasetT[k, j, i]
                        data_array.InsertNextValue(valueT)

            image_data.GetCellData().SetScalars(data_array)
            image_data.SetSpacing(1.0, 1.0, 1.0)
            image_data.SetOrigin(0.0, 0.0, 0.0)

            writer = vtk.vtkXMLImageDataWriter()
            writer.SetCompressorTypeToNone()
            writer.SetDataModeToAscii()
            writer.SetFileName(output_file)
            writer.SetInputData(image_data)
            writer.Write()

# Specify the HDF5 file path
file_path = "50.h5"

# Specify the steps to convert
steps = range(0, 1001, 100)

# Convert each step to VTI file
for step in steps:
    output_file = "C_Step{}.vti".format(step)
    convert_phi_to_vti(file_path, step, output_file)

#for step in steps:
    output_file = "T_Step{}.vti".format(step)
    convert_T_to_vti(file_path, step, output_file)
