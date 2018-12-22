using Makie, AbstractPlotting

using AbstractPlotting: RefValue, absrect
scene = scatter(rand(10), rand(10))

const key = Mouse.left

# Create an initially hidden rectangle

function select_rectangle(scene)

    waspressed = false # Why cant't I have waspressed be a normal Julia variable???
    rect = FRect() # Why can't I have `rect` be a normal Julia variable???


    rect_vis = lines!(
        scene,
        rect,
        linestyle = :dot,
        linewidth = 0.1,
        color = (:black, 0.4),
        visible = false,
        raw = true
    ).plots[end] # Why do I have to do .plots[end] ???


    # What does this line do? why do I need `camera`?
    # what is the `drag` argument after `do` ?
    dragged_rect = on(camera(scene), events(scene).mousedrag, key) do drag, key


        if ispressed(scene, key) && is_mouseinside(scene)
            mp = mouseposition(scene)
            if drag == Mouse.down
                waspressed = true
                rect_vis[:visible] = true # start displaying
                rect = FRect(mp, 0.0, 0.0)

                # What does this line do???
                rect_vis[1] = rect

            elseif drag == Mouse.pressed
                mini = minimum(rect[])
                rect = FRect(mini, mp - mini)

                # What does this line do???
                rect_vis[1] = rect
            end

        else
            if drag == Mouse.up && waspressed # User has selected the rectangle
                waspressed[] = false
                r = absrect(rect)
                w, h = widths(r)


                if w > 0.0 && h > 0.0 # Ensure that the rectangle has non0 size.
                    rect = r
                end
            end
            # always hide if not the right key is pressed
            rect_vis[:visible] = false # make the plotted rectangle invisible
        end

        return rect
    end
    return rect
end

dragged_rect = select_rectangle(scene)

on(dragged_rect) do
    println("rect = $(dragged_rect.origin) + $(dragged_rect.widths)")
end
