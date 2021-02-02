import os
import glob

configfile: 'config.yaml'


audio_file_paths = glob.glob("source_files/*.mp3")
audio_files = []

for audio_file_path in audio_file_paths:
    audio_files.append(os.path.basename(audio_file_path).split(".")[0])

def get_chunks(wildcards):
    checkpoint_output = checkpoints.chunk_audio.get(**wildcards).output[0]
    return expand("API_output/{audio_file}/chunk_{i}.json",
           audio_file=wildcards.audio_file,
           i=glob_wildcards(os.path.join(checkpoint_output, "chunk_{i}.wav")).i)


def get_chunks_wav(wildcards):
    checkpoint_output = checkpoints.chunk_audio.get(**wildcards).output[0]
    return expand("chunked_files/{audio_file}/chunk_{i}.wav",
           audio_file=wildcards.audio_file,
           i=glob_wildcards(os.path.join(checkpoint_output, "chunk_{i}.wav")).i)


rule all:
    input:
        #expand("API_output/{audio_file}/finished.txt", audio_file=audio_files)
        expand("output/{audio_file}.mp4", audio_file=audio_files),
        expand("output/{audio_file}.srt", audio_file=audio_files)

rule finish_diarization:
    input: expand("API_output/{audio_file}/finished.txt", audio_file=audio_files)



rule convert:
    input:
        "source_files/{audio_file}.mp3"
    output:
        temp("converted_files/{audio_file}.wav")
    shell:
        """ffmpeg -i {input} {output}"""

rule subsample:
    input:
        "converted_files/{audio_file}.wav"
    output:
        temp("converted_files/{audio_file}_16k.wav")
    shell:
        """sox {input} -b 16 {output} rate 16k"""

checkpoint chunk_audio:
    input:
        "converted_files/{audio_file}_16k.wav"
    output:
        directory("chunked_files/{audio_file}")
    shell:
        """mkdir {output}; ffmpeg -i {input} -f segment -segment_time 200 -c copy {output}/chunk_%03d.wav"""


rule send_chunks:
    input:
        "chunked_files/{audio_file}/chunk_{i}.wav"
    output:
        'API_output/{audio_file}/chunk_{i}.json'
    params:
        api_key = config['api_key']
    script:
        'scripts/api_conn.py'

rule finish:
    input:
        get_chunks
    output:
        "API_output/{audio_file}/finished.txt"
    shell:
        """touch {output}"""

rule convert_to_video:
    input:
        "source_files/{audio_file}.mp3"
    output:
        "output/{audio_file}.mp4"
    params:
        background="images/black_background.png"
    shell:
        """ffmpeg -loop 1 -i {params.background} -i {input} -c:v libx264 -tune stillimage -shortest {output}"""

rule flag_subtitles:
    input:
        json = get_chunks,
        wav = get_chunks_wav
    output:
        "output/{audio_file}.srt"
    script:
        "scripts/text_to_speech.py"
