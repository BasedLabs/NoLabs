<template>
  <q-drawer show-if-above :width="drawerWidth">
    <q-list>
      <q-item-label class="text-bold text-white" header>BioBuddy Chat</q-item-label>
      <q-separator />
      <q-item v-for="(messageGroup, index) in messages" :key="index">
        <q-item-section>
          <div class="q-pb-sm q-pt-sm text-bold">{{ displayName(messageGroup.role) }}</div>
          <div v-for="(message, msgIndex) in messageGroup.messages" :key="msgIndex">
            <div v-if="message.type === 'function'" class="function-message">
              <q-item-label class="text-h7 q-mb-sm">
                <div v-for="(functionCall, fcIndex) in message.content" :key="`function-${fcIndex}`">
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
                  <div class="markdown-content" v-html="renderMarkdown(message.content)"></div>
                  <q-btn flat v-if="messageGroup.role === 'user'" size="sm" icon="edit"
                    @click="startEditing(index)"></q-btn>
                  <q-btn flat v-if="isLastUserMessageWithoutResponse(index) && !sending" label="Regenerate"
                    icon-right="autorenew" @click="sendMessage"></q-btn>
                </div>
                <div v-else>
                  <q-input v-model="editMessage" dense filled></q-input>
                  <q-btn flat label="Save and Submit" @click="saveEdit(index)"></q-btn>
                  <q-btn flat label="Cancel" @click="cancelEdit"></q-btn>
                </div>
              </q-item-label>
            </div>
          </div>
          <q-spinner v-if="messageGroup.role === 'biobuddy' && messageGroup.awaitingActions" size="20px" />
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
import { FunctionCall, type FunctionParam, Message, GetAvailableFunctionCallsResponse } from 'src/refinedApi/client';
import { editMessageApi, loadConversationApi, saveMessageApi, sendQueryApi, getToolsApi } from 'src/features/biobuddy/api';

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
      messages: [] as { id: string, role: string, messages: { type: string, content: any }[], awaitingActions?: boolean }[],
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
      tools: null as GetAvailableFunctionCallsResponse | null,
      awaitingResponse: false,
      pendingActions: [] as string[],
    };
  },
  methods: {
    connectWebSocket() {
      this.socket = new WebSocket('ws://127.0.0.1:5738/ws');
      this.socket.onopen = () => {
        console.log('Connected to WebSocket server');
      };
      this.socket.onmessage = async (event) => {
        const message = JSON.parse(event.data);
        if (message.reply_type === 'stream') {
          this.currentMessageBuffer += message.content;
          const actionMatches = this.currentMessageBuffer.match(/<ACTION>(.*?)<END_ACTION>/g);
          if (actionMatches) {
            actionMatches.forEach((match) => {
              const actionText = match.replace(/<\/?ACTION>/g, '');
              this.pendingActions.push(actionText);
              this.currentMessageBuffer = this.currentMessageBuffer.replace(/<ACTION>/g, '');
              this.currentMessageBuffer = this.currentMessageBuffer.replace(/<END_ACTION>/g, '');
            });
          }
          if (this.awaitingResponse) {
            const lastMessageIndex = this.messages.length - 1;
            if (
              this.messages[lastMessageIndex].role === 'biobuddy' &&
              this.messages[lastMessageIndex].messages[this.messages[lastMessageIndex].messages.length - 1].type !== 'function'
            ) {
              this.messages[lastMessageIndex].messages[this.messages[lastMessageIndex].messages.length - 1].content = this.currentMessageBuffer;
            } else if (this.messages[lastMessageIndex].role === 'biobuddy') {
              this.messages[lastMessageIndex].messages.push({
                type: 'text',
                content: this.currentMessageBuffer,
              });
            } else if (this.messages[lastMessageIndex].role === 'user') {
              const newBioBuddyMessage = {
                id: new Date().getTime().toString(), // Generate a temporary ID
                role: 'biobuddy',
                messages: [
                  {
                    type: 'text',
                    content: this.currentMessageBuffer,
                  },
                ],
              };
              this.messages.push(newBioBuddyMessage);
            }
          } else {
            this.messages.push({
              id: new Date().getTime().toString(),
              role: 'biobuddy',
              messages: [{
                type: 'text',
                content: this.currentMessageBuffer,
              }]
            });
            this.awaitingResponse = true;
          }
        } else if (message.reply_type === 'final' || message.content.includes('<STOP>')) {
          this.awaitingResponse = false;
          const actionMatches = this.currentMessageBuffer.match(/<ACTION>(.*?)<END_ACTION>/g);
          if (actionMatches) {
            actionMatches.forEach((match) => {
              const actionText = match.replace(/<\/?ACTION>/g, '');
              this.pendingActions.push(actionText);
              this.currentMessageBuffer = this.currentMessageBuffer.replace(/<ACTION>/g, '');
              this.currentMessageBuffer = this.currentMessageBuffer.replace(/<END_ACTION>/g, '');
            });
          }
          const lastMessageIndex = this.messages.length - 1;
          if (this.messages[lastMessageIndex].role === 'biobuddy') {
            this.messages[lastMessageIndex].awaitingActions = true;
            await this.processPendingActions(this.messages[lastMessageIndex]);
          }
          this.saveMessage({
            id: new Date().getTime().toString(), // Generate a temporary ID
            role: 'assistant',
            messages: [{
              type: 'text',
              content: this.currentMessageBuffer,
            }],
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
    async processPendingActions(messageGroup) {
      for (const action of this.pendingActions) {
        const response = await sendQueryApi(
          action
        );
        messageGroup.messages.push({
          type: 'function',
          content: response.biobuddy_response.message,
        });
      }
      messageGroup.awaitingActions = false;
      this.pendingActions = [];
    },
    async loadConversation() {
      if (!this.experimentId) return;

      const response = await loadConversationApi(this.experimentId);
      const fetchedMessages = response.messages;

      const groupedMessages = [] as { id: string, role: string, messages: { type: string, content: any }[], awaitingActions?: boolean }[];
      let currentGroup = null as { id: string, role: string, messages: { type: string, content: any }[], awaitingActions?: boolean } | null;

      fetchedMessages.forEach((msg, index) => {
        if (msg.role === 'biobuddy') {
          if (!currentGroup) {
            currentGroup = {
              id: msg.id,
              role: 'biobuddy',
              messages: []
            };
          }
          currentGroup?.messages.push({
            type: msg.type,
            content: msg.message.content,
          });

          // If it's the last message, push the current group
          if (index === fetchedMessages.length - 1) {
            groupedMessages.push(currentGroup);
          }
        } else {
          // Push the current group if it exists
          if (currentGroup) {
            groupedMessages.push(currentGroup);
            currentGroup = null;
          }
          groupedMessages.push({
            id: msg.id,
            role: msg.role,
            messages: [{
              type: msg.type,
              content: msg.message.content,
            }]
          });
        }
      });

      this.messages = groupedMessages;
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
        messages: [
          {
            type: 'text',
            content: this.newMessage,
          },
        ],
      };

      this.messages.push(userMessage);
      await this.saveMessage(userMessage);

      const messagePayload = {
        experiment_id: this.experimentId,
        message_content: this.newMessage,
        previous_messages: this.messages.flatMap(msg =>
          msg.messages.map(m => ({
            role: msg.role,
            content: typeof m.content === 'string' ? m.content : JSON.stringify(m.content)
          }))
        ),
        tools: this.tools?.function_calls,
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
      for (const msg of message.messages) {
        if (msg.type === 'text') {
          await saveMessageApi(
            this.experimentId,
            msg.content,
            message.role,
          );
        } else if (msg.type === 'function') {
          // Assuming function call messages also need to be saved, adjust as necessary
          await saveMessageApi(
            this.experimentId,
            JSON.stringify(msg.content),
            message.role,
          );
        }
      }
    },
    isLastUserMessageWithoutResponse(index: number) {
      return index === this.messages.length - 1 && this.messages[index].role === 'user';
    },
    startEditing(index: number) {
      this.editIndex = index;
      this.editMessage = this.messages[index].messages.map(m => m.content).join(' ');
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
      this.messages[index].messages = [{
        type: 'text',
        content: this.editMessage,
      }];
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
    this.tools = await getToolsApi();
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
