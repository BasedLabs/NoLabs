<template>
  <q-drawer show-if-above :width="drawerWidth">
    <q-list>
      <q-item-label class="text-bold text-white" header>BioBuddy Chat</q-item-label>
      <q-separator />
      <q-item v-for="(message, index) in messages" :key="index">
        <q-item-section>
          <div v-if="message.type === 'function'" class="function-message">
            <q-item-label class="text-h7 q-mb-sm">
              <div class="q-pb-sm q-pt-sm text-bold">{{ displayName(message.role) }}</div>
              <div v-for="(functionCall, fcIndex) in message.message" :key="`function-${fcIndex}`">
                <p>
                  <q-icon name="check_circle" color="purple"></q-icon>
                  <span class="text-h7 text-purple q-ml-sm">{{ functionCall.function_name }}</span>
                </p>
                <ul>
                  Params:
                  <li v-for="(param, pIndex) in functionCall.arguments" :key="`param-${fcIndex}-${pIndex}`">
                    {{ param.name }}: {{ param.value }}
                  </li>
                </ul>
              </div>
            </q-item-label>
          </div>
          <div v-else>
            <q-item-label class="text-h7 q-mb-sm">
              <div v-if="editIndex !== index">
                <div class="q-pb-sm q-pt-sm text-bold">{{ displayName(message.role) }}</div>
                <div class="markdown-content" v-html="renderMarkdown(message.message.content)"></div>
                <q-btn flat v-if="message.role === 'user'" size="sm" icon="edit" @click="startEditing(index)"></q-btn>
                <q-btn flat v-if="isLastUserMessageWithoutResponse(index) && !sending" label="Regenerate" icon-right="autorenew" @click="sendQuery(message.message.content)"></q-btn>
              </div>
              <div v-else>
                <q-input v-model="editMessage" dense filled></q-input>
                <q-btn flat label="Save and Submit" @click="saveEdit(index)"></q-btn>
                <q-btn flat label="Cancel" @click="cancelEdit"></q-btn>
              </div>
            </q-item-label>
          </div>
        </q-item-section>
      </q-item>
    </q-list>
    <div class="q-pa-md">
      <q-input v-model="newMessage" label="Type a message..." dense filled @keyup.enter="sendMessage">
        <template v-slot:append>
          <q-btn icon="stop" flat @click="stopGeneration" :disable="!awaitingResponse || sending">Stop</q-btn>
          <q-btn icon="send" flat @click="sendMessage" :disable="sending">
            <q-spinner size="20px" v-if="sending"></q-spinner>
          </q-btn>
        </template>
      </q-input>
    </div>

    <div class="drawer-resize-handle" @mousedown="startResizing"></div>
  </q-drawer>
</template>

<script lang="ts">
import { Notify } from 'quasar';
import { defineComponent } from 'vue';
import MarkdownIt from 'markdown-it';
import { useBioBuddyStore } from './storage';
import { FunctionCall, type FunctionParam, Message } from 'src/refinedApi/client';
import { editMessageApi, loadConversationApi, saveMessageApi, sendQueryApi } from 'src/features/biobuddy/api';

export interface FunctionMapping {
  name: string;
  function: (parameters: any) => void;
}

export default defineComponent({
  name: 'BioBuddyChat',
  props: {
    experimentId: {
      type: String,
      required: true,
    },
  },
  data() {
    return {
      drawer: true,
      messages: [] as Message[],
      newMessage: '',
      sending: false,
      functionMappings: [] as FunctionMapping[],
      drawerWidth: 500, // Initial width of the drawer
      isResizing: false,
      initialMouseX: 0,
      editIndex: -1,
      editMessage: '',
      socket: null as WebSocket | null,
      currentMessageBuffer: '' as string,
      awaitingResponse: false,
    };
  },
  methods: {
    connectWebSocket() {
      this.socket = new WebSocket('ws://127.0.0.1:5738/ws');
      this.socket.onopen = () => {
        console.log('Connected to WebSocket server');
      };
      this.socket.onmessage = (event) => {
        const message = JSON.parse(event.data);
        if (message.reply_type === 'stream') {
          this.currentMessageBuffer += message.content;
          if (this.awaitingResponse) {
            const lastMessageIndex = this.messages.length - 1;
            if (this.messages[lastMessageIndex].role === 'biobuddy') {
              this.messages[lastMessageIndex].message.content = this.currentMessageBuffer;
            } else {
              this.messages.push({
                id: new Date().getTime().toString(), // Generate a temporary ID
                role: 'biobuddy',
                type: 'text',
                message: {
                  content: this.currentMessageBuffer,
                },
              });
              this.awaitingResponse = true;
            }
          } else {
            this.messages.push({
              id: new Date().getTime().toString(), // Generate a temporary ID
              role: 'biobuddy',
              type: 'text',
              message: {
                content: this.currentMessageBuffer,
              },
            });
            this.awaitingResponse = true;
          }
        } else if (message.reply_type === 'final' || message.content.includes('<STOP>')) {
          this.awaitingResponse = false;
          this.saveMessage({
            id: new Date().getTime().toString(), // Generate a temporary ID
            role: 'assistant',
            type: 'text',
            message: {
              content: this.currentMessageBuffer,
            },
          });
          this.currentMessageBuffer = '';
        }
      };
      this.socket.onclose = () => {
        console.log('WebSocket connection closed');
      };
      this.socket.onerror = (error) => {
        console.error('WebSocket error:', error);
      };
    },
    async loadConversation() {
      if (!this.experimentId) return;
      const response = await loadConversationApi(this.experimentId);
      this.messages = response.messages;
    },
    renderMarkdown(text: string) {
      const md = new MarkdownIt();
      let result = md.render(text);
      result = result.replace(
        /<h1>/g,
        '<h1 style="font-size: 1.5em; margin-top: 0.5em; margin-bottom: 0.5em; line-height: 1.2; font-weight: bold;">'
      );
      result = result.replace(/<h2>/g, '<h2 style="font-size: 1.3em; margin-top: 0.3em; margin-bottom: 0.3em;">');
      result = result.replace(/<h3>/g, '<h3 style="font-size: 1.1em; margin-top: 0.3em; margin-bottom: 0.3em;">');
      return result;
    },
    async sendMessage() {
      if (!this.experimentId || !this.newMessage.trim()) return;

      const userMessage = {
        id: new Date().getTime().toString(), // Generate a temporary ID
        role: 'user',
        type: 'text',
        message: {
          content: this.newMessage,
        },
      };

      this.messages.push(userMessage);
      await this.saveMessage(userMessage);

      const messagePayload = {
        experiment_id: this.experimentId,
        message_content: this.newMessage,
        previous_messages: this.messages.map(msg => ({ role: msg.role, content: msg.message.content })),
        tools: this.functionMappings,
        job_id: null, // Optional job ID
      };

      this.socket?.send(JSON.stringify(messagePayload));
      this.newMessage = '';
      this.awaitingResponse = true;
      this.currentMessageBuffer = ''; // Reset the buffer for new messages
    },
    async stopGeneration() {
      const stopPayload = {
        action: 'stop',
        experiment_id: this.experimentId,
      };
      this.socket?.send(JSON.stringify(stopPayload));
    },
    async saveMessage(message: any) {
      await saveMessageApi(
       this.experimentId,
       message.message.content,
       message.role,
      );
    },
    isLastUserMessageWithoutResponse(index: number) {
      return index === this.messages.length - 1 && this.messages[index].role === 'user';
    },
    startEditing(index: number) {
      this.editIndex = index;
      this.editMessage = this.messages[index].message.content;
    },
    cancelEdit() {
      this.editIndex = -1;
      this.editMessage = '';
    },
    async saveEdit(index: number) {
      if (!this.experimentId || !this.editMessage.trim()) return;
      const messageId = this.messages[index].id; // Assuming each message has a unique ID

      const editMessage = {
        experiment_id: this.experimentId,
        message_id: messageId,
        message_content: this.editMessage,
      };

      this.socket?.send(JSON.stringify(editMessage));
      this.messages[index].message.content = this.editMessage;
      this.messages = this.messages.slice(0, index + 1); // Remove messages after the edited one
      this.editIndex = -1;
      this.editMessage = '';
    },
    displayName(role: string) {
      return role === 'user' ? 'You' : 'Biobuddy';
    },
    startResizing(event: MouseEvent) {
      this.isResizing = true;
      this.initialMouseX = event.clientX;
      window.addEventListener('mousemove', this.resizeDrawer);
      window.addEventListener('mouseup', this.stopResizing);
    },
    resizeDrawer(event: MouseEvent) {
      if (!this.isResizing) return;
      const deltaX = event.clientX - this.initialMouseX;
      this.drawerWidth += deltaX;
      this.initialMouseX = event.clientX; // Update initial X to the current position
    },
    stopResizing() {
      this.isResizing = false;
      window.removeEventListener('mousemove', this.resizeDrawer);
      window.removeEventListener('mouseup', this.stopResizing);
    },
  },
  async mounted() {
    await this.loadConversation();
    this.connectWebSocket();
    const bioBuddyStore = useBioBuddyStore();
    this.functionMappings = [
      { name: 'query_rcsb_pdb_by_id', function: bioBuddyStore.invokeQueryRcsbPdbEventHandlers },
      { name: 'query_rcsb_pdb_by_protein_names', function: bioBuddyStore.invokeQueryRcsbPdbEventHandlers },
      { name: 'query_chembl', function: bioBuddyStore.invokeQueryChemblEventHandlers },
      { name: 'query_chembl_by_condition', function: bioBuddyStore.invokeQueryChemblEventHandlers },
    ];
  },
});
</script>

<style scoped>
.drawer-resize-handle {
  cursor: ew-resize;
  position: absolute;
  top: 0;
  right: 0; /* Adjust based on your layout, ensuring it's reachable for resizing */
  width: 20px;
  height: 100%;
}
</style>
