"""
Example usage of the discs_in_context package.
"""

from discs_in_context import plotcloud as pc
import matplotlib.pyplot as plt

# Example 1: Using a preset region (Taurus)
plot1 = pc(region='taurus', coord_system='icrs')
fig1, ax1 = plot1.plot(
    save_path='taurus_icrs.pdf',
    plot_pms=True,plot_discs=True
)

# Example 4: Different region (Lupus)
plot4 = pc(region='lupus', coord_system='galactic')
fig4, ax4 = plot4.plot(
    save_path='lupus_galactic.pdf',
    plot_pms=True
)

