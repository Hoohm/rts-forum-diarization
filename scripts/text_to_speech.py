import speech_recognition as sr
import json
from datetime import timedelta
import srt


def generate_srt(
    start_index,
    start_time,
    end_time,
    transcript,
    max_words_per_section=12,
    words_per_seconds=3,
):
    srt_parts = []
    words = transcript.split()
    if len(words) < 12:
        srt_section = srt.Subtitle(
            index=start_index, start=start_time, end=end_time, content=transcript
        )
        srt_parts.append(srt_section)
    else:

        return []


def get_chunk_subtitles(
    speaker_map_path, audio_in_path, start_index=1, offset=0, lag_correction=-0.001
):

    with open(speaker_map_path, "r") as in_json:
        data = json.load(in_json)
        results = data["channels"]["0"]
        n_speakers = int(results["summary"]["speaker-map"]["speaker_count"])
        raw_data = results["transitions"]["speaker-map"]

    audio_data = sr.AudioFile(audio_in_path)
    chunk_srt = []
    with audio_data as source:
        for entry in raw_data:

            offset_timedelta = timedelta(milliseconds=offset)

            start_time = (
                timedelta(milliseconds=(entry["timestamp_start"])) + offset_timedelta
            )
            end_time = (
                timedelta(milliseconds=(entry["timestamp_end"])) + offset_timedelta
            )

            speaker = entry["result"]
            duration = entry["timestamp_end"] / 1000 - entry["timestamp_start"] / 1000

            audio = r.record(source, duration=duration)

            # We only deal with speakers,
            if speaker != "no_speech":
                transcript = r.recognize_google(audio, language="fr-FR", show_all=True)
                # We take the first entry
                if len(transcript) != 0:
                    text = "{}: {}".format(
                        speaker, transcript["alternative"][0]["transcript"]
                    )

                    # generate_srt(
                    #     start_index,
                    #     start_time,
                    #     end_time,
                    #     transcript["alternative"][0]["transcript"],
                    # )
                    srt_part = srt.Subtitle(
                        index=start_index, start=start_time, end=end_time, content=text
                    )
                    chunk_srt.append(srt_part)
                    start_index += 1
            offset += lag_correction
    return start_index, entry["timestamp_end"], chunk_srt


r = sr.Recognizer()
start_index = 1
offset = 0
with open(snakemake.output[0], "w") as out_srt:
    for file_in, audio_in_path in zip(snakemake.input["json"], snakemake.input["wav"]):
        start_index, offset, chunk_srt = get_chunk_subtitles(
            speaker_map_path=file_in,
            audio_in_path=audio_in_path,
            start_index=start_index,
            offset=offset,
        )
        for chunk in chunk_srt:
            out_srt.write(chunk.to_srt())
