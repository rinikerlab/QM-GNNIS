import torch

class TinyModel(torch.nn.Module):

    def __init__(self, natoms):
        super(TinyModel, self).__init__()
        self.linear1 = torch.nn.Linear(3 * natoms, 200)
        self.activation = torch.nn.ReLU()

    def forward(self, coordinates):
        x = coordinates.flatten()
        x = self.linear1(x)
        x = self.activation(x)
        return torch.sum(x) * (10 ** -7)