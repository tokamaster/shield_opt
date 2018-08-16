import openmc
import openmc.model
import glob
import os

def shield(rad, bool, outer, m_batches, m_error, neutrons):

    #MATERIALS#

    min = 1
    max = outer

    R1 = rad[0]
    R2 = rad[1]

    mats = openmc.Materials()

    tungsten = openmc.Material(name='Tungsten')
    tungsten.set_density('g/cm3', 19.0)
    tungsten.add_element('W', 1.0)
    mats.append(tungsten)

    water = openmc.Material(name='Water')
    water.set_density('g/cm3', 1.0)
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    mats.append(water)

    #GEOMETRY#

    sphere1 = openmc.Sphere(R=min)
    sphere2 = openmc.Sphere(R=R1)
    sphere3 = openmc.Sphere(R=R2)
    sphere4 = openmc.Sphere(R=max, boundary_type='vacuum')

    vac1 = -sphere1
    mat1 = +sphere1 & -sphere2
    mat2 = +sphere2 & -sphere3
    mat3 = +sphere3 & -sphere4
    vac2 = +sphere4

    vacuum1 = openmc.Cell(region=vac1)
    first = openmc.Cell(region=mat1)
    first.fill = tungsten
    second = openmc.Cell(region=mat2)
    second.fill = water
    third = openmc.Cell(region=mat3)
    third.fill = tungsten
    vacuum2 = openmc.Cell(region=vac2)

    root = openmc.Universe(cells=(vacuum1, first, second, third,vacuum2))
    geom = openmc.Geometry(root)

    #SETTINGS#

    batch = 2
    max_batches = m_batches
    inactive = 1
    particles = neutrons

    source = openmc.Source()
    source.space = openmc.stats.Point((0,0,0))
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14e6], [1])

    sett = openmc.Settings()
    sett.batches = batch
    sett.trigger_active = True
    sett.trigger_max_batches = m_batches
    sett.inactive = inactive
    sett.particles = particles
    sett.output = {'tallies': False}
    sett.run_mode = 'fixed source'
    sett.source = source

    #TALLIES#

    tallies = openmc.Tallies()

    filter = openmc.SurfaceFilter(sphere4)

    tally = openmc.Tally(name='leakage')
    tally.scores = ['current']
    tally.filters = [filter]
    trigger = openmc.Trigger('rel_err', m_error/100)
    tally.triggers = [trigger]
    tallies.append(tally)

    model = openmc.model.Model(geom, mats, sett, tallies)

    #RUN#

    model.run(output=bool)

    #PLOT#

    #plots = openmc.Plots()

    #plot = openmc.Plot()
    #plot.basis = 'xz'
    #plot.origin = (0, 0, 0)
    #plot.width = (200, 200)
    #plot.pixels = (400, 400)
    #plot.color_by = 'material'
    #plot.colors = {tungsten: 'black', water: 'blue'}
    #plots.append(plot)

    #plots.export_to_xml()

    for i in reversed(range(batch,max_batches+1)):
        filename = 'statepoint.' +str(i).zfill(len(str(max_batches)))+ '.h5'
        if os.path.isfile(filename):
            sp = openmc.StatePoint(filename)
            break

    #print('file:',filename)
    leakage = sp.get_tally(name='leakage')
    leakage_mean = leakage.mean[0][0][0]
    leakage_error = leakage.std_dev[0][0][0]
    #print(leakage_mean)
    #print(leakage_error)

    dir = '/home/emiralle/shield_opt'
    for zippath in glob.iglob(os.path.join(dir, '*.h5')):
        os.remove(zippath)

    return leakage_mean*100, leakage_error*100
