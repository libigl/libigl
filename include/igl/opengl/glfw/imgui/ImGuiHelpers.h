// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Jérémie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_OPENGL_GLFW_IMGUI_IMGUIHELPERS_H
#define IGL_OPENGL_GLFW_IMGUI_IMGUIHELPERS_H

////////////////////////////////////////////////////////////////////////////////
#include <imgui/imgui.h>
#include <vector>
#include <string>
////////////////////////////////////////////////////////////////////////////////

// Extend ImGui by populating its namespace directly
//
// Code snippets taken from there:
// https://eliasdaler.github.io/using-imgui-with-sfml-pt2/
namespace ImGui
{

static auto vector_getter = [](void* vec, int idx, const char** out_text)
{
	auto& vector = *static_cast<std::vector<std::string>*>(vec);
	if (idx < 0 || idx >= static_cast<int>(vector.size())) { return false; }
	*out_text = vector.at(idx).c_str();
	return true;
};

inline bool Combo(const char* label, int* idx, std::vector<std::string>& values)
{
	if (values.empty()) { return false; }
	return Combo(label, idx, vector_getter,
		static_cast<void*>(&values), values.size());
}

inline bool ListBox(const char* label, int* idx, std::vector<std::string>& values)
{
	if (values.empty()) { return false; }
	return ListBox(label, idx, vector_getter,
		static_cast<void*>(&values), values.size());
}

} // namespace ImGui

#endif // IGL_OPENGL_GLFW_IMGUI_IMGUIHELPERS_H
