import http.client
import time
import json
import sys


def fetch_data(request_type, body_string, data, headers, conn):
    """Sends an http request and waits for a response. Will try and send 5 times if not successful

    args:
        request_type(str): Can be POST or GET
        body_string(str): The request itself
        data(str): ??
        headers(dict): A dict containing more information for the API.
        conn(HTTPSConnection): The connection to the API.
    """
    tries = 0
    response = ""
    for _ in range(5):
        time.sleep(1)
        conn.request(
            request_type,
            body_string,
            data,
            headers,
        )
        res = conn.getresponse()
        response = res.read()
        if (
            response
            == b"<html>\r\n<head><title>502 Bad Gateway</title></head>\r\n<body>\r\n<center><h1>502 Bad Gateway</h1></center>\r\n</body>\r\n</html>\r\n"
            or response is None
            or response == b""
        ):
            print(response)
            print("Attempt failed. trying again")
            time.sleep(8)
            tries += 1
        if tries == 2:
            sys.exit("Tried {} times. Exiting".format(tries))
    return response


# REPLACE with your key
API_KEY = snakemake.params.api_key
conn = http.client.HTTPSConnection("api.oto.ai")
# The input here is the wav file
with open(snakemake.input[0], "rb") as f:
    data = f.read()

post_headers = {"content-type": "audio/wav", "x-api-key": API_KEY}
post_body_string = "/file-processing/jobs?models=speaker-map&output_period=4096&include_summary=true&include_transitions=true&volume_threshold=0.001"

sent_request_response = fetch_data(
    request_type="POST",
    body_string=post_body_string,
    data=data,
    headers=post_headers,
    conn=conn,
)

job_id = json.loads(sent_request_response)["id"]
print("Job sent and id is: {}".format(job_id))
get_headers = {"x-api-key": API_KEY}

state = "new"
result = {}
check_body_string = "/file-processing/jobs/{}".format(job_id)
check_progress_response = fetch_data(
    request_type="GET",
    body_string=check_body_string,
    data="",
    headers=get_headers,
    conn=conn,
)
while not (state == "done" or state == "error"):
    time.sleep(1)
    check_progress_response = fetch_data(
        request_type="GET",
        body_string=check_body_string,
        data="",
        headers=get_headers,
        conn=conn,
    )
    result = json.loads(check_progress_response)
    state = result["state"]
    print("Job configuration:", result["config"])
    print("Job state:", state)
    if state == "error":
        print("There was an error: ", result["error_description"])
        sys.exit()
    if state == "done":
        time.sleep(1)
        fetch_body_string = "/file-processing/jobs/{}/results".format(job_id)
        final_data = fetch_data(
            request_type="GET",
            data="",
            body_string=fetch_body_string,
            headers=get_headers,
            conn=conn,
        )
        result = json.loads(final_data)["result"]
    time.sleep(30)
with open(snakemake.output[0], "w") as outfile:
    json.dump(result, outfile)