"""
Example usage of the discs_in_context package.
"""

from discs_in_context import plotcloud as pc
import matplotlib.pyplot as plt

# Example 1: Using a preset region (Taurus)
plot1 = pc(region='taurus', coord_system='icrs')
plot1.plot(
    save_path='taurus_icrs.pdf',
    plot_discs=True,plot_pms=True,interactive=True,dustmap='planck'
)

# Example 2: Different region (Lupus)
plot2 = pc(region='lupus', coord_system='galactic')
fig2, ax2 = plot2.plot(
    save_path='lupus_galactic.pdf',plot_discs=True,
    interactive=False
)

# Example 3: around a specific object
plot3 = pc(object='AB Aur', coord_system='icrs')
fig3, ax3 = plot3.plot(
    save_path='ab_aur_icrs.pdf',plot_discs=True,plot_pms=True,interactive=False,dustmap='planck'
)

# Example 4: all sky view
plot4 = pc(region='allsky',coord_system='galactic')
fig4, ax4 = plot4.plot(
    save_path='all_sky_galactic.pdf',plot_discs=True,interactive=True,dustmap='planck',show_regions=True
)

