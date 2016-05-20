% for enum in enums:
py::enum_<\
% for n in enum['namespaces']:
${n}::\
% endfor
${enum['name']}>(m, "${enum['name']}")
% for c in enum['constants']:
    .value("${c}", \
% for n in enum['namespaces']:
${n}::\
% endfor
${c})
% endfor
    .export_values();
% endfor


% for func in functions:
m.def("${func['name']}", []
(
  % for p in func['parameters'][:-1]:
  ${p['type']} ${p['name']},
  % endfor
  ${func['parameters'][-1]['type']} ${func['parameters'][-1]['name']}
)
{
  return \
% for n in func['namespaces']:
${n}::\
% endfor
${func['name']}(\
% for p in func['parameters'][:-1]:
${p['name']}, \
% endfor
${func['parameters'][-1]['name']});
}, __doc_\
% for n in func['namespaces']:
${n}_\
% endfor
${func['name']},
% for p in func['parameters'][:-1]:
py::arg("${p['name']}"), \
% endfor
py::arg("${func['parameters'][-1]['name']}"));

% endfor
