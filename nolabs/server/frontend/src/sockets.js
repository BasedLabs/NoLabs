import { io } from "socket.io-client";
import { baseUrl } from './storageConstants';

const socket = io(baseUrl);

export default socket;