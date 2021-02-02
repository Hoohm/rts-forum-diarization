# Purpose

This pipeline allows to use the [OTO API](https://docs.oto.ai/api/)

Here are the main steps:
1. Take audio files located in the `source_files`.
2. convert and chunk into the right size for the API.
3. Send the chunks and retrieve the speaker transition timings.
4. Use the speech to text google API to get the text from the audio.
5. Generate an `.srt` subtitle file
6. Convert the original mp3 into an mp4 with black background
7. Enjoy your subbed speaker detected dialog!

![dag](images/dag.png?raw = true "dag")

# Configuration

The `template/config.yaml` holds all the important configuration fields.

* api_key: Your API key for OTO.ai

# Aknowledgements

I would like to thank personnaly Nicolas Peroni for giving me access to the speaker diarization BETA.
Starting this project is an old dream of mine and this API is getting me one step closer.
I would also like to thank all the *[OTO.ai]*(https://www.oto.ai/) team for working on this tool and helping me out with the few bumps on the road!