#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuizmoMenu.h>
#include <imgui/imgui.h>
#include <igl/readSTL.h>
#include <igl/material_colors.h>
#include <iostream>
#include "tutorial_shared_path.h"

class SlicingPlugin : public igl::opengl::glfw::imgui::ImGuiMenu {

	igl::opengl::ViewerData data;

	virtual void init(igl::opengl::glfw::Viewer *_viewer) {
		using namespace igl;
		ImGuiMenu::init(_viewer);

		// Load slicing plane into viewer (a square mesh)
		const Eigen::MatrixXd V = (Eigen::MatrixXd(4, 3) <<
			-0.5, -0.5, 0.0,
			-0.5,  0.5, 0.0,
			 0.5,  0.5, 0.0,
			 0.5, -0.5, 0.0).finished();
		const Eigen::MatrixXi F = (Eigen::MatrixXi(2, 3) <<
			0, 2, 1,
			0, 3, 2).finished();

		data.set_mesh(V, F);
		data.set_face_based(true);
		data.set_colors(Eigen::RowVector4d(224, 86, 253, 128)/255.0);
		data.show_lines = false;
	}

	virtual bool pre_draw() override {
		ImGuiMenu::pre_draw();
		ImGuizmo::BeginFrame();
		ImGui::Begin("ImGuizmo Tools");
		return false;
	}

	virtual bool post_draw() override {
		viewer->core().draw(data);
		ImGuiMenu::post_draw();
		return false;
	}

	virtual void draw_custom_window() override {

		static Eigen::Matrix4f matrix = Eigen::Matrix4f::Identity();
		Eigen::Matrix4f proj = viewer->core().proj;
		Eigen::Affine3f rescale = Eigen::Scaling(0.5f * viewer->core().camera_base_zoom)
			* Eigen::Translation3f(viewer->core().camera_base_translation);
		Eigen::Affine3f view = Eigen::Affine3f(viewer->core().view*(1./viewer->core().camera_zoom)) 
                            * rescale.inverse();
		
    // Pass Viewer transformation matrices to ImGuizmo
    ImGuizmo::EditTransform(view.matrix().data(), proj.data(), matrix.data());

    // Transform the slicing plane according to 
    // ImGuizmo tool manipulations in the viewer
		Eigen::Affine3f model(matrix);
		model = rescale.inverse() * model;

		Eigen::MatrixXd V = (Eigen::MatrixXd(4, 3) <<
			-0.5, -0.5, 0.0,
			-0.5,  0.5, 0.0,
			 0.5,  0.5, 0.0,
			 0.5, -0.5, 0.0).finished();
		V = (V.rowwise().homogeneous() 
        * model.matrix().cast<double>().transpose()).rowwise().hnormalized();

		const Eigen::MatrixXi F = (Eigen::MatrixXi(2, 3) <<
			0, 2, 1,
			0, 3, 2).finished();

		if (data.V.rows() == V.rows())
			data.set_vertices(V);
	}
};

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd N;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  std::string filename(TUTORIAL_SHARED_PATH "/cow.off");
  viewer.load_mesh_from_file(filename);

  // Custom menu
  SlicingPlugin menu;
  viewer.plugins.push_back(&menu);

  viewer.resize(1024, 1024);
  viewer.launch();
}