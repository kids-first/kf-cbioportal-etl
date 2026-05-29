"""Helper functions for downloading files from URLs with retry logic and parallelization.

Overall function created from co-pilot prompt, refined and debugged by M. Brown
"""

import logging
import random
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep

import requests
from requests.adapters import HTTPAdapter

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

adapter = HTTPAdapter(pool_connections=20, pool_maxsize=20)

def _request_with_retries(
        session: requests.Session,
        url: str,
        headers: dict[str, str] | None = None,
        stream: bool = True,
        retries: int = 5,
        delay: int = 3
        ) -> requests.Response| None:
    """Perform a GET request with retries and exponential backoff.

    Args:
        session (requests.Session): The requests session to use.
        url (str): The URL to request.
        headers (dict[str, str], optional): Headers to include in the request. Defaults to None.
        stream (bool, optional): Whether to stream the response. Defaults to True.
        retries (int, optional): Number of retry attempts on failure. Defaults to 5.
        delay (int, optional): Initial delay between retries in seconds. Defaults to 3.

    Returns:
        requests.Response | None: The response object or None if all retries fail.

    """
    logger = logging.getLogger(__name__)

    for attempt in range(retries):
        try:
            response =  session.get(url, headers=headers, stream=stream, timeout=30)
            response.raise_for_status()

            if attempt > 0:
                logger.info("Successfully downloaded %s on attempt %d", url, attempt + 1)
            return response
        except (requests.RequestException, requests.Timeout, requests.HTTPError) as e:  # noqa: PERF203
            logger.warning("Error downloading %s", url)
            if attempt < retries - 1:
                logger.warning("Retrying in %d seconds...", delay)
                sleep(delay + random.uniform(0, 1))
                delay *= 2  # Exponential backoff
            else:
                msg = f"Failed to download {url} after {retries} attempts because of error: {e}"
                raise Exception(msg) from e
    return None


def download_range(session: requests.Session,
                   url: str,
                   start: int,
                   end: int,
                   retries=5,
                   delay=3) -> tuple[int, bytearray]:
    """Download a specific byte range from a URL with retry logic.

    Args:
        session (requests.Session): The requests session to use.
        url (str): The URL to download from.
        start (int): The starting byte position.
        end (int): The ending byte position.
        retries (int): Number of retry attempts on failure.
        delay (int): Initial delay between retries in seconds.

    Returns:
        tuple[int, bytes]: A tuple containing the starting byte position and the downloaded data.

    """
    logger = logging.getLogger(__name__)
    headers = {"Range": f"bytes={start}-{end}"}
    with _request_with_retries(session, url, headers, True, retries, delay) as response:
        # Validate range support
        if response.status_code != 206:
            v_err = f"Expected 206 Partial Content, got {response.status_code}"
            raise ValueError(v_err)

        content_range = response.headers.get("Content-Range")
        logger.debug("Content-Range: %s", content_range)
        expected_prefix = f"bytes {start}-{end}/"

        if not content_range or not content_range.startswith(expected_prefix):
            v_err = f"Invalid Content-Range: {content_range}, expected {expected_prefix}"
            raise ValueError(v_err)
        data = bytearray()
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                data.extend(chunk)
    # Integrity check
    logger.debug("Downloaded chunk size: %d, expected: %d", len(data), end - start + 1)
    expected_size = end - start + 1
    if len(data) != expected_size:
        v_err = f"Incomplete chunk {start}-{end}: got {len(data)}, expected {expected_size}"
        raise ValueError(v_err)

    return start, data


def parallel_download(
        url: str,
        output_path: str,
        total_size: int,
        num_workers: int = 8,
        chunk_size: int = 32 * 1024 * 1024
        ) -> None:
    """Download a file in parallel using multiple threads.

    Args:
        url (str): The URL to download from.
        output_path (str): The path to save the downloaded file.
        total_size (int): The total size of the file in bytes.
        num_workers (int): The number of parallel threads to use.
        chunk_size (int): The size of each chunk to download.

    """
    logger = logging.getLogger(__name__)
    ranges = [
        (start, min(start + chunk_size - 1, total_size - 1))
        for start in range(0, total_size, chunk_size)
        ]
    logger.debug("Download ranges: %s", ranges)
    with open(output_path, "wb") as f:
        f.truncate(total_size)  # preallocate

    session = requests.Session()
    session.mount("https://", adapter)

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(download_range, session, url, start, end)
            for start, end in ranges
        ]

        with open(output_path, "r+b") as f:
            for future in as_completed(futures):
                start, data = future.result()
                f.seek(start)
                f.write(data)

    session.close()


def small_download(url: str, output_path: str, retries: int = 5, delay: int = 3) -> None:
    """Download a small file with retry logic.

    Args:
        url (str): The URL to download from.
        output_path (str): The path to save the downloaded file.
        retries (int): Number of retry attempts on failure.
        delay (int): Initial delay between retries in seconds.

    """
    session = requests.Session()
    session.mount("https://", adapter)

    response = _request_with_retries(session, url, None, False, retries, delay)
    with open(output_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)

    session.close()
