#include <imguizmo/ImGuizmo.h>
#include "menu_extra.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/material_colors.h>
#include <imgui/imgui.h>
#include <iostream>


class SlicingPlugin : public igl::opengl::glfw::imgui::ImGuiMenu {

	igl::opengl::ViewerData data;

	virtual void init(igl::opengl::glfw::Viewer *_viewer) {
		using namespace igl;
		ImGuiMenu::init(_viewer);

		// Inline mesh of a square
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
		return false;
	}

	virtual bool post_draw() override {
		viewer->core().draw(data);
		ImGuiMenu::post_draw();
		return false;
	}

	virtual void draw_custom_window() override {
		static Eigen::Matrix4f matrix = Eigen::Matrix4f::Identity();
		Eigen::Affine3f rescale = Eigen::Scaling(0.5f * viewer->core().camera_base_zoom)
			* Eigen::Translation3f(viewer->core().camera_base_translation);
		// Eigen::Affine3f view = Eigen::Affine3f(viewer->core().view) * rescale.inverse();
		Eigen::Matrix4f proj = viewer->core().proj;

		float camera_zoom = viewer->core().camera_zoom;
		Eigen::Affine3f view = Eigen::Affine3f(viewer->core().view * (1./camera_zoom)) * rescale.inverse();

		ImGuizmo::EditTransform(view.matrix().data(), proj.data(), matrix.data());

		Eigen::Affine3f model(matrix);
		model = rescale.inverse() * model;
		// std::cout << viewer->core().camera_base_zoom << std::endl;
		// std::cout << viewer->core().camera_base_translation.transpose() << std::endl;

		Eigen::MatrixXd V = (Eigen::MatrixXd(4, 3) <<
			-0.5, -0.5, 0.0,
			-0.5,  0.5, 0.0,
			 0.5,  0.5, 0.0,
			 0.5, -0.5, 0.0).finished();
		V = (V.rowwise().homogeneous() * model.matrix().cast<double>().transpose()).rowwise().hnormalized();

		const Eigen::MatrixXi F = (Eigen::MatrixXi(2, 3) <<
			0, 2, 1,
			0, 3, 2).finished();

		// data.clear();
		if (data.V.rows() == V.rows())
			data.set_vertices(V);
		// data.set_face_based(true);
		// data.set_colors(Eigen::RowVector4d(149, 165, 166, 128)/255.0);
		// data.show_lines = false;
	}

};
