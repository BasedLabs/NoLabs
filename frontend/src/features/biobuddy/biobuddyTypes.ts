// workflowTypes.ts

// Define the type for each connection in a component
export interface BiobuddyComponentConnection {
  source_component_id: string;
  source_output: string;
  target_input: string;
}

// Define the type for each component
export interface BiobuddyWorkflowComponent {
  id: string;
  name: string;
  description: string;
  x: number;
  y: number;
  connections: BiobuddyComponentConnection[];
}

// Define the main data type
export interface BiobuddyWorkflowAdjustmentData {
  workflow_components: BiobuddyWorkflowComponent[];
}
