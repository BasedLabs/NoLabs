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
          <div v-if="messageGroup.role === 'biobuddy' && messageGroup.awaitingActions" class="q-mt-sm q-mb-sm">
            <q-spinner size="20px" />
            <span class="text-h7 text-purple q-ml-sm">{{ currentAction }}</span>
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
import { v4 as uuidv4 } from 'uuid';
import { useBioBuddyStore } from './storage';
import { type FunctionParam, FunctionCall_Output, Message, GetAvailableFunctionCallsResponse } from 'src/refinedApi/client';
import { editMessageApi, loadConversationApi, saveMessageApi, sendQueryApi, getToolsApi, saveFunctionCallApi } from 'src/features/biobuddy/api';
import { MeasurementFlags } from 'ngl/dist/declarations/component/structure-component';

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
      messages: [] as { id: string, role: string, messages: { id: string, type: string, content: any }[], awaitingActions?: boolean }[],
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
      currentAction: ''
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
              const actionText = match.replace(/<\/?ACTION>/g, '').replace(/<\/?END_ACTION>/g, '');;
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
                id: uuidv4(),
                type: 'text',
                content: this.currentMessageBuffer,
              });
            } else if (this.messages[lastMessageIndex].role === 'user') {
              const newBioBuddyMessage = {
                id: uuidv4(), // Generate a temporary ID
                role: 'biobuddy',
                messages: [
                  {
                    id: uuidv4(),
                    type: 'text',
                    content: this.currentMessageBuffer,
                  },
                ],
              };
              this.messages.push(newBioBuddyMessage);
            }
          } else {
            const newMessage = {
              id: uuidv4(),
              type: 'text',
              content: this.currentMessageBuffer,
            }
            this.messages.push({
              id: uuidv4(),
              role: 'biobuddy',
              messages: [newMessage]
            });
            this.saveMessage(newMessage);
            this.awaitingResponse = true;
          }
        } else if (message.reply_type === 'final' || message.content.includes('<STOP>')) {
          this.awaitingResponse = false;
          const actionMatches = this.currentMessageBuffer.match(/<ACTION>(.*?)<END_ACTION>/g);
          if (actionMatches) {
            actionMatches.forEach((match) => {
              const actionText = match.replace(/<\/?ACTION>/g, '').replace(/<\/?END_ACTION>/g, '');
              this.pendingActions.push(actionText);
              this.currentMessageBuffer = this.currentMessageBuffer.replace(/<ACTION>/g, '');
              this.currentMessageBuffer = this.currentMessageBuffer.replace(/<END_ACTION>/g, '');
            });
          }
          const lastMessageIndex = this.messages.length - 1;
          this.saveMessage({
            id: uuidv4(),
            type: 'text',
            role: 'biobuddy',
            content: this.currentMessageBuffer,
          });
          if (this.messages[lastMessageIndex].role === 'biobuddy') {
            this.messages[lastMessageIndex].awaitingActions = true;
            await this.processPendingActions(this.messages[lastMessageIndex]);
          }
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
        this.currentAction = action;
        const response = await sendQueryApi(
          this.experimentId,
          action
        );
        const biobuddyMessage = {
          id: uuidv4(),
          type: 'function',
          role: 'biobuddy',
          content: response.biobuddy_response.message,
        }
        messageGroup.messages.push(biobuddyMessage);
        this.saveMessage(biobuddyMessage)
      }
      messageGroup.awaitingActions = false;
      this.pendingActions = [];
    },
    async loadConversation() {
      if (!this.experimentId) return;

      const response = await loadConversationApi(this.experimentId);
      const fetchedMessages = response.messages;

      const groupedMessages = [] as { id: string, role: string, messages: { id: string, type: string, content: any }[], awaitingActions?: boolean }[];
      let currentGroup = null as { id: string, role: string, messages: { id: string, type: string, content: any }[], awaitingActions?: boolean } | null;

      fetchedMessages.forEach((msg, index) => {
        if (msg.role === 'biobuddy' && msg.type === 'text') {
          if (!currentGroup) {
            currentGroup = {
              id: msg.id,
              role: 'biobuddy',
              messages: []
            };
          }
          currentGroup?.messages.push({
            id: msg.message.id,
            type: msg.type,
            content: msg.message.content,
          });
          // If it's the last message, push the current group
          if (index === fetchedMessages.length - 1) {
            groupedMessages.push(currentGroup);
          }
        } else if (msg.role === 'biobuddy' && msg.type === 'function') {
          if (!currentGroup) {
            currentGroup = {
              id: msg.id,
              role: 'biobuddy',
              messages: []
            };
          }
          currentGroup?.messages.push({
            id: msg.id,
            type: msg.type,
            content: msg.message,
          });
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
              id: msg.message.id,
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
        id: uuidv4(), // Generate a temporary ID
        role: 'user',
        messages: [
          {
            id: uuidv4(),
            type: 'text',
            role: 'user',
            content: this.newMessage,
          },
        ],
      };

      this.messages.push(userMessage);
      for (const msg of userMessage.messages) {
        await this.saveMessage(msg);
      }

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
    async saveMessage(msg: any) {
      if (msg.type === 'text') {
        await saveMessageApi(
          this.experimentId,
          msg.id,
          msg.content,
          msg.role,
        );
      } else if (msg.type === 'function') {
        // Assuming function call messages also need to be saved, adjust as necessary
        await saveFunctionCallApi(
          this.experimentId,
          msg.id,
          msg.content[0],
          msg.role
        );
        await this.invokeFunction(msg.content[0]);
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

      this.messages[index].messages = [{
        id: messageId,
        type: 'text',
        content: this.editMessage,
      }];
      this.messages = this.messages.slice(0, index + 1);

      const messagePayload = {
        experiment_id: this.experimentId,
        message_content: this.editMessage,
        previous_messages: this.messages.flatMap(msg =>
          msg.messages.map(m => ({
            role: msg.role,
            content: typeof m.content === 'string' ? m.content : JSON.stringify(m.content)
          }))
        ),
        tools: this.tools?.function_calls,
        job_id: null, // Optional job ID
      };

      this.socket?.send(JSON.stringify(messagePayload)); // Remove messages after the edited one
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
    async invokeFunction(functionCall: FunctionCall_Output) {
      const mapping = this.functionMappings?.find(m => m.name === functionCall.function_name);
      if (mapping && typeof mapping.function === 'function') {
        await mapping.function(functionCall.data);
      }
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
