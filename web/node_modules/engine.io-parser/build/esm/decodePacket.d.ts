import { Packet, BinaryType, RawData } from "engine.io-parser/build/esm/commons.js";
declare const decodePacket: (encodedPacket: RawData, binaryType?: BinaryType) => Packet;
export default decodePacket;
