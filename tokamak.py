import openmc
import openmc.model
import glob
import os
from numpy import sqrt

def shield(rad, major, coil, bool, outer, m_batches, m_error, neutrons):

    #MATERIALS#

    min = coil
    max = outer+coil
    major_rad = major + outer

    R1 = rad[0]
    R2 = rad[1]

    mats = openmc.Materials()

    tungsten = openmc.Material(name='Tungsten carbide')
    tungsten.set_density('g/cm3', 15.63)
    tungsten.add_element('W', 1.0)
    tungsten.add_element('C', 1.0)
    mats.append(tungsten)

    water = openmc.Material(name='Water')
    water.set_density('g/cm3', 1.0)
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    mats.append(water)

    copper = openmc.Material(name='Copper')
    copper.set_density('g/cm3', 8.5)
    copper.add_element('Cu', 1.0)
    mats.append(copper)

    eurofer = openmc.Material(name='EUROFER97')
    eurofer.set_density('g/cm3', 7.75)
    eurofer.add_element('Fe', 89.067, percent_type='wo')
    eurofer.add_element('C', 0.11, percent_type='wo')
    eurofer.add_element('Mn', 0.4, percent_type='wo')
    eurofer.add_element('Cr', 9.0, percent_type='wo')
    eurofer.add_element('Ta', 0.12, percent_type='wo')
    eurofer.add_element('W', 1.1, percent_type='wo')
    eurofer.add_element('N', 0.003, percent_type='wo')
    eurofer.add_element('V', 0.2, percent_type='wo')
    mats.append(eurofer)

    #GEOMETRY#

    cylinder = openmc.ZCylinder(R=min)

    shield1 = openmc.ZCylinder(R=R1)
    shield2 = openmc.ZCylinder(R=R2)
    shield3 = openmc.ZCylinder(R=max)
    top = major_rad-(max-min)
    shield_sph1 = openmc.Sphere(R=top)
    shield_sph2 = openmc.Sphere(R=major_rad-(max-min)+(max-R2))
    shield_sph3 = openmc.Sphere(R=major_rad-(max-min)+(max-R2)+(R2-R1))

    pit = sqrt(top**2+top**2)
    max_shield = sqrt(top**2+pit**2)+1

    vessel_in = openmc.Sphere(R=major_rad)
    vessel_out = openmc.Sphere(R=max_shield, boundary_type='vacuum')

    magnet = -cylinder & -vessel_in

    mat1 = -shield1 & +cylinder & -shield_sph3  #work on this
    mat2 = -shield2 & +shield1 & -shield_sph2 #work on this
    mat3 = -shield3 & +shield2 & -shield_sph1 #work on this

    plasma = +shield3 & -shield_sph1

    a = +shield_sph1 & -shield_sph2 & +shield2 #work on this
    b = +shield_sph2 & -shield_sph3 & +shield1 #work on this
    c = +shield_sph3 & -vessel_in & +cylinder #work on this

    vessel = +vessel_in & -vessel_out

    plasma_cell = openmc.Cell(region=plasma) #1

    cop_mag = openmc.Cell(region=magnet) #2
    cop_mag.fill = copper

    tung1 = openmc.Cell(region=mat1) #3
    tung1.fill = tungsten
    wat = openmc.Cell(region=mat2) #4
    wat.fill = water
    tung2 = openmc.Cell(region=mat3) #5
    tung2.fill = tungsten

    ste_ves = openmc.Cell(region=vessel) #6
    ste_ves.fill = eurofer

    tung_sph1 = openmc.Cell(region=a) #7
    tung_sph1.fill = tungsten
    wat_sph = openmc.Cell(region=b) #8
    wat_sph.fill = water
    tung_sph2 = openmc.Cell(region=c) #9
    tung_sph2.fill = tungsten

    root = openmc.Universe(cells=(plasma_cell, cop_mag, tung1, wat, tung2, ste_ves, tung_sph1, wat_sph, tung_sph2))
    geom = openmc.Geometry(root)
    root.plot(width=(600.0, 600.0), basis='xz')

    #SETTINGS#

    batch = 2
    max_batches = m_batches
    inactive = 0
    particles = neutrons

    source = openmc.Source()
    source.space = openmc.stats.Box((-top,-top,-top),(top,top,top))
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14e6], [1])

    sett = openmc.Settings()
    sett.batches = batch
    sett.trigger_active = True
    sett.trigger_max_batches = m_batches
    sett.inactive = inactive
    sett.particles = particles
    sett.output = {'tallies': True}
    sett.run_mode = 'fixed source'
    sett.source = source

    #PLOT#

    plots = openmc.Plots()

    plot = openmc.Plot()
    #plot.type = 'voxel'
    plot.basis = 'xz'
    plot.origin = (0, 0, 0)
    plot.width = (600, 600)
    plot.pixels = (1000, 1000)
    plot.color_by = 'material'
    plot.colors = {tungsten: 'black', water: 'blue', eurofer: 'grey', copper: 'brown'}
    plots.append(plot)

    plots.export_to_xml()

    #TALLIES#

    tallies = openmc.Tallies()

    filter = openmc.CellFilter(cop_mag)

    tally = openmc.Tally(name='total')
    tally.scores = ['total']
    tally.filters = [filter]
    trigger = openmc.Trigger('rel_err', m_error/100)
    tally.triggers = [trigger]
    tallies.append(tally)

    model = openmc.model.Model(geom, mats, sett, tallies)

    #RUN#

    model.run(output=bool)

    for i in reversed(range(batch,max_batches+1)):
        filename = 'statepoint.' +str(i).zfill(len(str(max_batches)))+ '.h5'
        if os.path.isfile(filename):
            sp = openmc.StatePoint(filename)
            break

    #print('file:',filename)
    leakage = sp.get_tally(name='total')
    leakage_mean = leakage.mean[0][0][0]
    leakage_error = leakage.std_dev[0][0][0]
    #print(leakage_mean)
    #print(leakage_error)

    dir = '/home/emiralle/shield_git'
    for zippath in glob.iglob(os.path.join(dir, '*.h5')):
        os.remove(zippath)

    return leakage_mean, leakage_error
