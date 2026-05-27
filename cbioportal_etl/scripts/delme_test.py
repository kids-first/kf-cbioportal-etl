from time import sleep

import sevenbridges as sbg
import sys
import requests

from concurrent.futures import ThreadPoolExecutor, as_completed



def download_range(url: str, start: int, end: int, retries=5, delay=3) -> tuple[int, bytes]:
    """Download a specific byte range from a URL with retry logic.

    Args:
        url (str): The URL to download from.
        start (int): The starting byte position.
        end (int): The ending byte position.
        retries (int): Number of retry attempts on failure.
        delay (int): Initial delay between retries in seconds.

    Returns:
        tuple[int, bytes]: A tuple containing the starting byte position and the downloaded data.

    """
    headers = {"Range": f"bytes={start}-{end}"}
    for attempt in range(retries):
        try:
            response = requests.get(url, headers=headers, stream=True, timeout=30)
            response.raise_for_status()

            data = bytearray()
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    data.extend(chunk)
            break
        except (requests.RequestException, requests.Timeout) as e:  # noqa: PERF203
            print(f"Error downloading range {start}-{end}: {e}", file=sys.stderr)
            if attempt < retries - 1:
                print(f"Retrying in {delay} seconds...", file=sys.stderr)
                sleep(delay)
                delay *= 2  # Exponential backoff
            else:
                msg = f"Failed to download range {start}-{end} after {retries} attempts"
                raise Exception(msg) from e
    return start, data


def parallel_download(url, output_path, total_size, num_workers=8):

    ranges = [
    (start, min(start + chunk_size - 1, total_size - 1))
    for start in range(0, total_size, chunk_size)
]

    for i in range(num_workers):
        start = i * chunk_size
        # Make sure last chunk reaches end of file
        end = (start + chunk_size - 1) if i < num_workers - 1 else total_size - 1
        ranges.append((start, end))

    with open(output_path, "wb") as f:
        f.truncate(total_size)  # preallocate

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(download_range, url, start, end)
            for start, end in ranges
        ]

        with open(output_path, "r+b") as f:
            for future in as_completed(futures):
                start, data = future.result()
                f.seek(start)
                f.write(data)

    print("Download complete!")


config: sbg.Config = sbg.Config(profile=sys.argv[1])
api = sbg.Api(config=config)

file_obj = api.files.get(sys.argv[2])

file_url = file_obj.download_info().url
chunk_size = 32 * 1024 * 1024
out = file_obj.name
total_size: int = file_obj.size
parallel_download(file_url, out, total_size, num_workers=16)

