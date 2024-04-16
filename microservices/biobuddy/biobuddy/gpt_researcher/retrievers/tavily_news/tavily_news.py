# Tavily API Retriever

# libraries
import os
from tavily import TavilyClient


class TavilyNews():
    """
    Tavily News API Retriever
    Retrieve news articles from the Tavily News API
    """
    def __init__(self, query):
        """
        Initializes the TavilySearch object
        Args:
            query:
        """
        self.query = query
        self.api_key = self.get_api_key()
        self.client = TavilyClient(self.api_key)

    def get_api_key(self):
        """
        Gets the Tavily API key
        Returns:

        """
        # Get the API key
        try:
            api_key = os.environ["TAVILY_API_KEY"]
        except:
            raise Exception("Tavily API key not found. Please set the TAVILY_API_KEY environment variable. "
                            "You can get a key at https://app.tavily.com")
        return api_key

    def search(self, max_results=7):
        """
        Searches the query
        Returns:

        """
        # Search the query
        results = self.client.search(self.query, search_depth="advanced", topic="news", max_results=max_results)
        # Return the results
        search_response = [{"href": obj["url"], "body": obj["content"]} for obj in results.get("results", [])]
        return search_response
