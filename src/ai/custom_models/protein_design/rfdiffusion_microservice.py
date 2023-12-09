from flask import Flask, request, jsonify

app = Flask(__name__)

@app.route('/rfdiffusion', methods=['POST'])
def rfdiffusion():
    try:
        # Assuming the incoming data is in JSON format
        data = request.get_json()

        # Accessing a key named 'message' from the JSON data
        message = data.get('message', 'No message provided.')

        # You can perform any processing with the received data here

        response = {'status': 'success', 'message': f'Message received: {message}'}

        return jsonify(response)

    except Exception as e:
        response = {'status': 'error', 'message': str(e)}
        return jsonify(response), 500  # Internal Server Error

if __name__ == '__main__':
    app.run(debug=True)