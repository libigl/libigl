#pragma once
////////////////////////////////////////////////////////////////////////////////
#include <imguizmo/ImGuizmo.h>
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>
#include <vector>
#include <Eigen/Dense>
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
