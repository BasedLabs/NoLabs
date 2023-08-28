import torch

class SimpleGOMultiLayerPerceptron(torch.nn.Module):
    def __init__(self, input_dim, num_classes):
        super(MultiLayerPerceptron, self).__init__()
        
        self.linear1 = torch.nn.Linear(input_dim, 1012)
        self.activation1 = torch.nn.ReLU()
        self.linear2 = torch.nn.Linear(1012, 840)
        self.activation2 = torch.nn.ReLU()
        self.linear3 = torch.nn.Linear(840, num_classes)
        
    def forward(self, x):
        x = self.linear1(x)
        x = self.activation1(x)
        x = self.linear2(x)
        x = self.activation2(x)
        x = self.linear3(x)
        return x

class SimpleSolubilityMultiLayerPerceptron(torch.nn.Module):
    # Current MLP is build on top of ESM-2 150M
    def __init__(self, input_size = 640, hidden_sizes = [320, 160, 80, 40], output_size = 1, dropout_rate = 0.2):
        super(MLPClassifier, self).__init__()
        layers = []
        prev_size = input_size
        for size in hidden_sizes:
            layers.append(nn.Linear(prev_size, size))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout_rate))  # Add dropout layer
            prev_size = size
        self.hidden_layers = nn.Sequential(*layers)
        self.output_layer = nn.Linear(prev_size, output_size)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.hidden_layers(x)
        x = self.output_layer(x)
        x = self.sigmoid(x)
        return x