import torch
import torch.nn as nn

class SimpleGOMultiLayerPerceptron(nn.Module):
    def __init__(self, input_dim, num_classes):
        super(SimpleGOMultiLayerPerceptron, self).__init__()
        self.linear1 = nn.Linear(input_dim, 500)
        self.activation1 = nn.ReLU()
        self.linear2 = nn.Linear(500, 300)
        self.activation2 = nn.ReLU()
        self.linear3 = nn.Linear(300, num_classes)
        
    def forward(self, x):
        x = self.linear1(x)
        x = self.activation1(x)
        x = self.linear2(x)
        x = self.activation2(x)
        x = self.linear3(x)
        return x

class SimpleSolubilityMultiLayerPerceptron(nn.Module):
    # Current MLP is build on top of ESM-2 150M
    def __init__(self, input_size = 640, hidden_sizes = [320, 160, 80, 40], output_size = 1, dropout_rate = 0.2):
        super(SimpleSolubilityMultiLayerPerceptron, self).__init__()
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