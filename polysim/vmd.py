#!/usr/bin/env python3
import subprocess

# Path to renderer
TACH = '/packages/apps/spack/18/opt/spack/gcc-12.1.0/vmd-1.9.4-khs/lib64/tachyon_LINUXAMD64'

def draw_step(trj, **kwargs):
    script = 'vmd-script.vmd'
    open(script, 'w').write(vmd_script(trj, **kwargs))
    cmd = ['vmd', '-dispdev', 'none', '-e', script]
    out = subprocess.run(cmd, capture_output=1, text=1)


def vmd_script(trj, **kwargs):
    s = [
    f'mol new {{{trj}}} type {{lammpstrj}} first 0 last -1 step 1 waitfor all',
    'display projection Orthographic',
    'display shadows on',
    'display ambientocclusion on',
    'display depthcue off',
    'color Display Background white',
    'axes location Off',
    'mol modstyle 0 0 VDW 0.2 12',
    'mol modcolor 0 0 Type',
    'mol material 0 0 AOChalky', 
    'color Type 1 lime' ]

    if 'step' in kwargs:
        s += [f'animate goto {kwargs["step"]}', 'display update']
    #,
    #'mol representation VDW 0.50 12',
    #'mol material AOChalky',
    #
    s += [
     'rotate x by 30',
     'rotate y by 30',
    f'render Tachyon vmdscene.dat "{TACH}" -aasamples 12 %s -format TARGA -o %s.tga',
    'quit']
    return '\n'.join(s)

