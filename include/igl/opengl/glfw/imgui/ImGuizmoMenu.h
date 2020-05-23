#pragma once
////////////////////////////////////////////////////////////////////////////////
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/material_colors.h>
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>
#include <imguizmo/ImGuizmo.h>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace ImGuizmo {
void EditTransform(const float *cameraView, float *cameraProjection, float* matrix)
{
	static ImGuizmo::OPERATION mCurrentGizmoOperation(ImGuizmo::ROTATE);
	static ImGuizmo::MODE mCurrentGizmoMode(ImGuizmo::WORLD);
	static bool useSnap = false;
	static float snap[3] = { 1.f, 1.f, 1.f };

	// if (ImGui::IsKeyPressed(90))
	// 	mCurrentGizmoOperation = ImGuizmo::TRANSLATE;
	// if (ImGui::IsKeyPressed(69))
	// 	mCurrentGizmoOperation = ImGuizmo::ROTATE;
	// if (ImGui::IsKeyPressed(82)) // r Key
	// 	mCurrentGizmoOperation = ImGuizmo::SCALE;
	if (ImGui::RadioButton("Translate", mCurrentGizmoOperation == ImGuizmo::TRANSLATE))
		mCurrentGizmoOperation = ImGuizmo::TRANSLATE;
	ImGui::SameLine();
	if (ImGui::RadioButton("Rotate", mCurrentGizmoOperation == ImGuizmo::ROTATE))
		mCurrentGizmoOperation = ImGuizmo::ROTATE;
	ImGui::SameLine();
	if (ImGui::RadioButton("Scale", mCurrentGizmoOperation == ImGuizmo::SCALE))
		mCurrentGizmoOperation = ImGuizmo::SCALE;
	float matrixTranslation[3], matrixRotation[3], matrixScale[3];
	ImGuizmo::DecomposeMatrixToComponents(matrix, matrixTranslation, matrixRotation, matrixScale);
	ImGui::InputFloat3("Tr", matrixTranslation, 3);
	ImGui::InputFloat3("Rt", matrixRotation, 3);
	ImGui::InputFloat3("Sc", matrixScale, 3);
	ImGuizmo::RecomposeMatrixFromComponents(matrixTranslation, matrixRotation, matrixScale, matrix);

	if (mCurrentGizmoOperation != ImGuizmo::SCALE)
	{
		if (ImGui::RadioButton("Local", mCurrentGizmoMode == ImGuizmo::LOCAL))
			mCurrentGizmoMode = ImGuizmo::LOCAL;
		ImGui::SameLine();
		if (ImGui::RadioButton("World", mCurrentGizmoMode == ImGuizmo::WORLD))
			mCurrentGizmoMode = ImGuizmo::WORLD;
	}
	// if (ImGui::IsKeyPressed(83))
	// 	useSnap = !useSnap;
	// ImGui::Checkbox("##Snap", &useSnap);
	// ImGui::SameLine();

	// switch (mCurrentGizmoOperation)
	// {
	// case ImGuizmo::TRANSLATE:
	// 	ImGui::InputFloat3("Snap", &snap[0]);
	// 	break;
	// case ImGuizmo::ROTATE:
	// 	ImGui::InputFloat("Angle Snap", &snap[0]);
	// 	break;
	// case ImGuizmo::SCALE:
	// 	ImGui::InputFloat("Scale Snap", &snap[0]);
	// 	break;
	// }

	ImGuiIO& io = ImGui::GetIO();
	ImGuizmo::SetRect(0, 0, io.DisplaySize.x, io.DisplaySize.y);
	ImGuizmo::Manipulate(cameraView, cameraProjection, mCurrentGizmoOperation, mCurrentGizmoMode, matrix, NULL, useSnap ? &snap[0] : NULL);
}
} // namespace ImGuizmo


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
		Eigen::Matrix4f proj = viewer->core().proj;

		float camera_zoom = viewer->core().camera_zoom;
		Eigen::Affine3f view = Eigen::Affine3f(viewer->core().view * (1./camera_zoom)) * rescale.inverse();

		ImGuizmo::EditTransform(view.matrix().data(), proj.data(), matrix.data());

		Eigen::Affine3f model(matrix);
		model = rescale.inverse() * model;

		Eigen::MatrixXd V = (Eigen::MatrixXd(4, 3) <<
			-0.5, -0.5, 0.0,
			-0.5,  0.5, 0.0,
			 0.5,  0.5, 0.0,
			 0.5, -0.5, 0.0).finished();
		V = (V.rowwise().homogeneous() * model.matrix().cast<double>().transpose()).rowwise().hnormalized();

		const Eigen::MatrixXi F = (Eigen::MatrixXi(2, 3) <<
			0, 2, 1,
			0, 3, 2).finished();

		if (data.V.rows() == V.rows())
			data.set_vertices(V);
	}

};
