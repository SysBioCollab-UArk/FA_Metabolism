import math


# Helper function for getting the optimal number of columns for a multi-plot figure
def get_fig_ncols(ndims):
    if not isinstance(ndims, int) or ndims <= 0:
        raise ValueError("'ndims' must be a positive integer")
    r1 = round(math.sqrt(ndims))  # round() returns an int
    r2 = math.ceil(math.sqrt(ndims))  # math.ceil() also returns an int
    while r1 * r2 >= ndims:
        r1 -= 1
        r2 += 1
    return min(r1 + 1, r2 - 1)  # the smaller of the two integers is the # of columns


# Function to hide extra plots at the end of a multi-plot figure and add x-axis labels and ticks to plots directly
# above the hidden plots
def hide_extra_plots(fig, axs_2D, n_plots):
    fig.canvas.draw()

    nrows, ncols = axs_2D.shape
    axs = axs_2D.flatten()

    xtick_fontsize = axs[n_plots - 1].get_xticklabels()[0].get_fontsize()
    xlabel_text = axs[n_plots - 1].get_xlabel()
    xlabel_fontsize = axs[n_plots - 1].xaxis.label.get_fontsize()

    visible_ticks_and_labels = [
        (val, label.get_text())
        for val, label in zip(axs[n_plots - 1].get_xticks(), axs[n_plots - 1].get_xticklabels())
        if axs[n_plots - 1].get_xlim()[0] <= val <= axs[n_plots - 1].get_xlim()[1]
    ]
    tick_vals, tick_labels = zip(*visible_ticks_and_labels)

    i = n_plots
    while i < nrows * ncols:
        axs[i].set_visible(False)

        # Get the data-to-display transform and the display-to-figure transform
        to_display = axs[i - ncols].transData.transform
        to_figure = fig.transFigure.inverted().transform

        # Determine y-position for tick labels just below axs[3]
        _, y_data = axs[i - ncols].get_ylim()
        y_display = axs[i - ncols].transAxes.transform((0, 0))[1]  # bottom of plot in display coords
        y_label = to_figure((0, y_display - 10))[1]  # shift down 25 pixels

        # Add tick labels manually under axs[3]
        for tick, label in zip(tick_vals, tick_labels):
            x_display = to_display((tick, 0))[0] - 6  # + 4 + 2 * (len(label) - 1) # add a small offset
            x_fig = to_figure((x_display, 0))[0]
            fig.text(x_fig, y_label, label, ha='center', va='top', fontsize=xtick_fontsize)

        # Add x-axis label centered under axs[3]
        x0, x1 = axs[i - ncols].get_position().x0, axs[i - ncols].get_position().x1
        xcenter = (x0 + x1) / 2
        fig.text(xcenter, y_label - 0.035, xlabel_text, ha='center', va='top', fontsize=xlabel_fontsize)

        # If this is the right-most column, add a 2nd y-axis to add space to the right, so last tick doesn't get cut off
        if i % ncols:
            twin_ax = axs[i - ncols].twinx()  # create a secondary y-axis on the right side
            twin_ax.set_ylabel(' ', fontsize=1, labelpad=12)  # add a dummy y-axis label
            twin_ax.tick_params(right=False, labelright=False)  # make tick marks and labels invisible

        # Move to next plot
        i += 1
