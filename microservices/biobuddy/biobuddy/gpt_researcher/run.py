from websocket import WebSocket

from .master.agent import GPTResearcher


async def get_report(query: str, websocket: WebSocket, report_type: str = 'research_report') -> str:
    """
    Get the report for the given query and report type:
    Args:
        query: The query to search for
        report_type: The type of report to generate (e.g. research_report)
    Output:
        report: The generated report in markdown format
    """
    researcher = GPTResearcher(query, report_type, config_path="../config.json", websocket=websocket)
    print("GPT researcher")
    report = await researcher.run()
    return report
