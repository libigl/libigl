#include <igl/opengl/glfw/Viewer.h>

bool load_rawfile(const std::string& rawfilename, const Eigen::RowVector3i& dims, Eigen::VectorXd& out, bool normalize=true) {
  const size_t num_bytes = dims[0]*dims[1]*dims[2];

  char* data = new char[num_bytes];
  std::ifstream rawfile(rawfilename, std::ifstream::binary);

  if (!rawfile.good()) {
    std::cerr << "ERROR: RawFile '" << rawfilename << "' does not exist." << std::endl;
    return false;
  }

  rawfile.read(data, num_bytes);
  if (!rawfile) {
    std::cerr << "ERROR: Only read " << rawfile.gcount() <<
                 " bytes from Raw File '" << rawfilename <<
                 "' but expected to read " << num_bytes <<
                 " bytes." << std::endl;
    return false;
  }
  rawfile.close();

  out.resize(num_bytes);
  for (int i = 0; i < num_bytes; i++) {
    out[i] = static_cast<double>(data[i]);
    if (normalize) {
      static_assert(sizeof(char) == sizeof(std::uint8_t), "Your system is fucked"); // This is dumb but why not
      out[i] /= 255.0;
    }
  }

  return true;
}

int main(int argc, char *argv[])
{
  // Inline mesh of a cube
  const Eigen::MatrixXd V= (Eigen::MatrixXd(8,3)<<
    0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    0.0,1.0,1.0,
    1.0,0.0,0.0,
    1.0,0.0,1.0,
    1.0,1.0,0.0,
    1.0,1.0,1.0).finished();
  const Eigen::MatrixXi F = (Eigen::MatrixXi(12,3)<<
    1,7,5,
    1,3,7,
    1,4,3,
    1,2,4,
    3,8,7,
    3,4,8,
    5,7,8,
    5,8,6,
    1,5,6,
    1,6,2,
    2,6,8,
    2,8,4).finished().array()-1;

  Eigen::VectorXd data;
  constexpr const char* file = "/home/francis/projects/fish_deformation/data/p-tapinosoma.raw";
  Eigen::RowVector3i dims = { 304, 234, 1719 };
  load_rawfile(file, dims, data, true);


  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.core.background_color = Eigen::RowVector4f(0., 0., 0., 1.);
  int volume_id = viewer.append_volume();
  viewer.volume().set_data(dims, data);

  viewer.data().set_mesh(V, F);
  viewer.data().clear();
  viewer.data().set_face_based(true);
  viewer.launch();
}
