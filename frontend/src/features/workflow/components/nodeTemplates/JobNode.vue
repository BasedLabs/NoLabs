<template>
  <q-card dark bordered :class="{'pulsating-border': isRunning}" v-if="!nodeData">
    <q-spinner color="primary" size="3em" />
  </q-card>
  <q-card dark bordered :class="{'pulsating-border': isRunning}" v-if="nodeData" style="max-width: 400px">
    <q-card-section v-if="nodeData?.data.error">
      <div class="row no-wrap items-center">
        <q-space />
        <q-btn v-if="nodeData?.data.error" color="red" class="q-pm-md">
          <q-icon left name="warning" />
          <div>{{ nodeData?.data.error }}</div>
        </q-btn>
      </div>
    </q-card-section>
    <q-card-section>
      <JobNodeContent :nodeId="nodeId" :name="nodeData?.name" :description="nodeData?.description"/>
    </q-card-section>
    <q-card-section>
      <NodeHandles :nodeId="nodeId" :inputs="nodeData?.data.inputs" :outputs="nodeData?.data.outputs" />
    </q-card-section>
    <q-card-actions align="around">
      <q-btn v-if="!isRunning" @click="startWorkflow" color="info" icon="play_arrow" label="Start" />
      <q-btn v-else>
        <q-spinner color="primary" size="20px" />
        Running...
      </q-btn>
      <q-btn @click="openDialogue" label="Extended view" />
      <q-btn @click="deleteNode" color="negative" icon="delete" label="Delete" />
    </q-card-actions>
  </q-card>
</template>


<script>
import JobNodeContent from './JobNodeContent.vue';
import NodeHandles from './NodeHandles.vue';
import { useWorkflowStore } from 'src/features/workflow/components/storage';

export default {
  name: 'JobNode',
  components: {
    JobNodeContent,
    NodeHandles
  },
  props: {
    nodeId: {
      type: String,
      required: true
    },
    onDeleteNode: Function,
    onOpenSettings: Function,
    onOpenDialog: Function
  },
  data() {
    return {
      nodeData: null
    };
  },
  computed: {
    isRunning() {
      const workflowStore = useWorkflowStore();
      return workflowStore.runningComponentIds.includes(this.nodeId);
    }
  },
  watch: {
    'useWorkflowStore().runningComponentIds': function(newVal, oldVal) {
      if (newVal !== oldVal) {
        this.updateRunningStatus();
      }
    }
  },
  mounted() {
    const workflowStore = useWorkflowStore();
    this.nodeData = workflowStore.getNodeById(this.nodeId);
    console.log(this.nodeData.name);
  },
  methods: {
    async startWorkflow() {
      const workflowStore = useWorkflowStore();
      const workflowId = workflowStore.workflowId;
      try {
        await startWorkflowComponent(workflowId, this.nodeId);
        this.isRunning = true;
        console.log(`Started workflow component with workflowId: ${workflowId} and nodeId: ${this.nodeId}`);
      } catch (error) {
        console.error('Failed to start workflow component', error);
      }
    },
    deleteNode() {
      if (this.onDeleteNode) {
        this.onDeleteNode(this.nodeId);
      }
    },
    openSettings() {
      if (this.onOpenSettings) {
        this.onOpenSettings(this.nodeId);
      }
    },
    openDialogue() {
      if (this.onOpenDialog) {
        this.onOpenDialog(this.nodeId);
      }
    },
    updateRunningStatus() {
      const workflowStore = useWorkflowStore();
      this.isRunning = workflowStore.runningComponentIds.includes(this.nodeId);
    }
  }
};
</script>

<style scoped>
.pulsating-border {
  border: 3px solid #42b983;
  animation: pulsate 1.5s ease-out infinite;
}

@keyframes pulsate {
  0% {
    box-shadow: 0 0 10px 2px rgba(66, 185, 131, 0.5);
  }
  50% {
    box-shadow: 0 0 10px 5px rgba(66, 185, 131, 1);
  }
  100% {
    box-shadow: 0 0 10px 2px rgba(66, 185, 131, 0.5);
  }
}
</style>
