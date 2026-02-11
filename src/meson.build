beng207_2_sources = [
  'main.c',
  'beng207_2-application.c',
  'beng207_2-window.c',
]

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : true)

beng207_2_deps = [
  dependency('gtk4'),
  dependency('libadwaita-1', version: '>= 1.4'),
  m_dep,
]

beng207_2_sources += gnome.compile_resources('beng207_2-resources',
  'beng207_2.gresource.xml',
  c_name: 'beng207_2'
)

executable('beng207_2', beng207_2_sources,
  dependencies: beng207_2_deps,
       install: true,
)