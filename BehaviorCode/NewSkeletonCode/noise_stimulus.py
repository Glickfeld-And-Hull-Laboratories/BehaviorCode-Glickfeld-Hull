import os
import time


# Absolute path to the directory where the noise images are stored
fname = getvar('folderName')
imagepath = ('/Users/hullglick/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/images/NoiseStimulusCapture/' + fname)

# Name of the subdirectory for the current session (date and time stamp)
sessiondir = '.' #time.strftime('%Y%m%d-%H%M%S')

# In order for MWorks to access the noise image directory, we need to create a
# symbolic link to it inside the current working directory
os.symlink(os.path.join(imagepath, sessiondir), 'images')


def capture_display_image(event):
    image_data = event.data

    # Ignore announcements that don't correspond to active captures or that
    # contain no data
    if not ('capture_display_image' in get_reverse_codec() and
            getvar('capture_display_image') and
            image_data):
        return

    # Create a per-trial subdirectory inside the session directory
    trialdir = 'trial_%d' % getvar('trial_number')
    trialpath = os.path.join(imagepath, sessiondir, trialdir)
    os.makedirs(trialpath, exist_ok=True)

    filename = os.path.join(trialpath, 'frame_%d.png' % getvar('image_number'))
    if os.path.isfile(filename):
        error('Image file already exists: ' + filename)
    else:
        with open(filename, 'wb') as fp:
            fp.write(event.data)
        message('Created image ' + filename)

    # Let the experiment know we're done
    setvar('capture_display_image', False)


# Remove any event callbacks left over from a previous experiment
unregister_event_callbacks()

# Register the capture callback
register_event_callback('#stimDisplayCapture', capture_display_image)
